function [cellInfos,adjMatrix] = findTree(cellList)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function cellStructure = findTree(cellList)
%@author:  Manuel Campos
%@date:    July 3, 2013
%@copyright 2012-2013 Yale University
%==========================================================================
%**********output********:
%cellInfos: Structure in which each entry of the structure correspond to a 
% cell during one cell cycle with the following information:
%       - oldID:        Cell id in the cellList, as in cellList.cellId
%       - frameArray:   array of frame numbers in which the cell is present
%       - ancestor:     Cell id of the mother cell (id as in the new
% structure)
%       - descendants:  Array of two cell ids of the two daughter cells (id
% as in the new structure)
%       - sister:       Cell id of the sister cell (id as in the new
% structure)
%adjMatrix:     Adjacency matrix for filiation relationship between cells
%**********Input********:
%cellList: cell structure that needs to be updated with extra fields.
%=========================================================================
% PURPOSE:  Transform a cellList into a structure containing one entry for
% each cell appearing throughout the time lapse. For each cell, all the
% frame in which the cell is detected are aggregated into one array. The
% same apply to cell length and area and a growth rate gr is fitted to the
% growth curve using the exponential model of growth L=Lo.exp(gr.time).
%The adjacency matrix for relatives allows for interrogating the cellInfos
%structure under a genealogy perspective
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% add extra fields to the cellList structure
% cellList = getExtraDataMultiThread(cellList);

% Clean cellList of bad cells (mesh=0)
tic
for f=1:length(cellList.meshData)
    frame=length(cellList.meshData)+1-f;
    cellNbFrame=length(cellList.meshData{frame});
    for c=1:cellNbFrame
        cell=cellNbFrame+1-c;
%         [frame, cell]
        goodCell=~isempty(cellList.meshData{frame}{cell}) && ...
            (isfield(cellList.meshData{frame}{cell},'mesh') && length(cellList.meshData{frame}{cell}.mesh)>1);
        if goodCell==0
            [frame, cell]
            cellList.meshData{frame}(cell)=[];
            cellList.cellId{frame}(cell)=[];
        end
    end
end
toc
%% Define presenceMat to vectorize the search for cell lifespans and genealogy
% Fill a matrix of 1 and 0 to specify the presence or absence of each cell
% in each frame

presenceMat=zeros(10000,length(cellList.meshData));
maxCellNd=zeros(length(cellList.meshData),1);
for frame = 1:length(cellList.meshData)
    [~, cellId] = CL_getFrame(frame,cellList);
    presenceMat(cellId,frame)=1;
    maxCellNd(frame)= max(cellId);
end
% Clean presenceMat (remove everything beyond max(maxCellNd))
% presenceMat=presenceMat(1:max(maxCellNd),:);
liveCells=find(sum(presenceMat,2)>0);

%% New cell IDs attribution in advance
% Attribute a number to each cell in advance so as to keep track of
% descendant and suster cell numbers as they are re-arrayed into the growth
% structure array. For any given cell undergoing division, a new number is
% associated to the cell after division (i.e. the mother cell id is
% different from both daughter cell ids).

newIDs=zeros(10000,3);count=1;
for ids = 1:size(liveCells,1)
    % determine the frame numbers in which cell ids is present
    indexArray=find(presenceMat(liveCells(ids),:)==1);
    % Get the cell info of this cell
    idCell=cellList.cellId{indexArray(end)}==liveCells(ids); % idCell is the id of the cell in the meshData{frame} cell array
    infoCell=cellList.meshData{indexArray(end)}{idCell};
    progeniture=[NaN,infoCell.descendants];
    for g=1:length(progeniture)
        newIDs(count,1)=count;
        newIDs(count,2)=liveCells(ids);
        newIDs(count,3)=progeniture(g);
        count=count+1;
    end
end
% Prepare the newIDs array for the next step
newIDs=newIDs(find(newIDs(:,1)~=0),:);
newIDs(end+1,:)=[NaN 0 NaN];

%% growth structure containing cell growth, metric and genealogy informations
% filled based on newIDs. The new cell ID correspond to the index of the
% old cell id in the mewIDs array. For cells undergoing divisions, the new
% ID corresponds to the index of (division event id + 1)th appearance of
% the old cell id in the newIDs array
% i.e. if cell one divides twice during the course of the experiment, the
% old id "1" will appear three times in newIDs array and the new ID for
% this cell before the first division will be the index of the first "1" in
% new IDs, the new ID for this cell after division 1 and before division 2
% will be the index of the second appearance of "1" in the newIDs array and
% so on.

for ids = 1:size(liveCells,1)
    % determine the frame numbers in which cell ids is present
    indexArray=find(presenceMat(liveCells(ids),:)==1);
    % Get the cell info of this cell
    idCell=cellList.cellId{indexArray(end)}==liveCells(ids);
    infoCell=cellList.meshData{indexArray(end)}{idCell};
    indices=find(newIDs(:,2)==liveCells(ids));
    idGen=find(newIDs(:,3)==liveCells(ids));

    cellEvents=[indexArray(1),infoCell.divisions,indexArray(end)];% birth, divisions, death
    infoCell.descendants(end+1)=0; % Facilitates the descendants search
    if isempty(infoCell.ancestors) % Facilitates ancestor search
        mother=NaN; 
    end
    
    for generations=1:length(indices)
        % Keep track of old cell ID
        cellInfos(indices(generations)).oldID=liveCells(ids);
        % determine the lifespan of the cell for this cell cycle
        frameSpan=logical((indexArray >= cellEvents(generations))-(indexArray >= cellEvents(generations+1)));
        cellInfos(indices(generations)).frameArray=indexArray(frameSpan);
        % Keep track of the ID of the mother cell
        cellInfos(indices(generations)).ancestor=mother;
        % Keep track of the ID of the sister cell
        cellInfos(indices(generations)).sister=idGen;
        % Keep track of the ID of the daughter cells
        index=[indices;NaN];indices=[indices;length(newIDs)];
        daughterOne = index(generations+1);
        daughterTwo = newIDs(find(newIDs(:,2)==infoCell.descendants(generations),1),1);
        cellInfos(indices(generations)).descendants=[daughterOne,daughterTwo];
        % Provide ancestor information for next generation
        mother=indices(generations);
    end
end

% Create tree of relatives. Ancestor in line, descendant in column
adjMatrix=zeros(length(cellInfos));
for i=1:length(cellInfos)
    if ~any(isnan(cellInfos(i).descendants))
        adjMatrix(i,cellInfos(i).descendants)=1;
%         if growth(i).descendants>726
%             i
%         end
    end
end
toc

%% Fill with data such as length, area, volume, compute growth rate
dt=0.12; % time interval in min for fitting purposes
for c=1:length(cellInfos)
    cellInfos(c).length=nan(1,length(cellInfos(c).frameArray));
    cellInfos(c).area=nan(1,length(cellInfos(c).frameArray));
    for frame=1:length(cellInfos(c).frameArray)
        cellStr=CL_getCell(cellInfos(c).oldID,cellInfos(c).frameArray(frame),cellList);
        cellInfos(c).length(frame)=cellStr.length;
        cellInfos(c).area(frame)=cellStr.area;
    end
    if length(cellInfos(c).frameArray)>9
        X=dt.*cellInfos(c).frameArray;
        Y=cellInfos(c).length;
%         options.startpoint=[growth(i).area(1) 0];
        [fitobject gof]=fit(X',Y','exp1');%,'options',fitStart);
        cellInfos(c).gr=fitobject.b;
        cellInfos(c).gof=gof;
    end
    c
end
toc
end
