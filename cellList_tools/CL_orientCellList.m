function cellList = CL_orientCellList(cellList,maxi,signal,pole,getExtraData,timelapse)
%function cellList = orientCellList(cellList,maxi,signal,pole,getExtraData,timelpase)
%microbeTracker vX.X.X
%@author: Géraldine Laloux; the reorient function 
%(reorient_parOfMicrTrck.m) and the kymograph function that are used 
%inside this function have been written by Oleksii Sliusarenko and must be 
%present in the path (reorient_partOfMicrTrck.m is in the same folder as
%this function).
%@date: March 21, 2013
%modified by:  Ahmad Paintdakhi - April 18, 2013.
%@copyright 2012-2013 Yale University
%==========================================================================
%*****************output******************:
%orientedList:      a new structure with similar organization as
%                   the input cellList (see below for changes applied in this new cellList)
%*****************input*******************:
%cellList:  the structure that needs to be reoriented (typically cellList)
%maxi:      optional - the maximum number of frames (from the time-lapse) that will be
%           taken into account for each cell. If [] is specified, the default value (the total
%           number of frames available) will be used.
%signal:    the signal (0,1,2,3,...) that will be used to orient the cells.
%           Write it as 'signal1', for example.
%pole:      0 for old pole, 1 for new pole. Determines whether the cell half where the most
%           signal is present should be considered as new pole or old
%           pole-proximal ("top" values or "bottom" values in the cell structure,
%           respectively). 
%==========================================================================
%Orients the cells based on the presence of more signal (due to the presence of a polar focus for
%example) in one half of the cell. For each cell where necessary, will flip the whole structure for
%that cell for each frame (until "maxi" or the last frame available
%depending on the value of maxi), so the "old pole", defined by the presence of more signal in one or the 
%other half of the cell (as specified by the 'pole' input), corresponds to
%the "bottom" part of the structure (and the new pole will be on top).
% This function can be useful if one wants to make sure that the old
% pole will always be placed at the bottom of single-cell kymographs from a time-lapse, 
% for example.
%**************************************************************************
%This function needs access to the library of microbeTracker.
%perform addpath(genpath('//aunt/common/microbeTrackerRepo')) if the 
%specified path and all its subfolders are not already in the matlab path.
%**************************************************************************
%==========================================================================   


% variables
cellToFlip = [];
totalCells = [];
if getExtraData == 1
	if ~isfield(cellList,'meshData'),  cellList = CL_makeCellListNew(cellList); end
	if isfield(cellList,'meshData'), cellList = CL_makeDouble(cellList); end
	for ii = 1:length(cellList.meshData)
		for jj = 1:length(cellList.meshData{ii})
			cellList.meshData{ii}{jj} = getextradata(cellList.meshData{ii}{jj});
		end
	end
end
if timelapse
    [~,cellIds] = CL_getFrame(1,cellList);
    for cell = cellIds
        cellStructure = CL_getCell(cell,1,cellList);
    % Gets last frame before division (determined by size of kymograph)          
        if isfield (cellStructure,signal)
           kymo = kymograph(cellList,cell,'signal0','nodisp');
           sizeKymo = size(kymo,2); 
    % Determines which frame will be the last one to be taken into account for this particular cell
    %(size of kymo or maxi, whichever is the smallest)           
            if sizeKymo < maxi
               lastFrame = sizeKymo;
            else
               lastFrame = maxi;
            end
    % Calculates fluorescence intensity for input signal in each half of the cell. Only cells with
    % signal information will be included in the new structure.
            if ~isempty (cellStructure.(signal))
               totalCells = [totalCells cell]; %#ok<AGROW>
               cellLength = length(cellStructure.lengthvector);
               halfCell = round(cellLength/2);
               halfCell1 = halfCell+1;
               bottomFluo = sum(cellStructure.(signal)(1:halfCell));
               upperFluo = sum(cellStructure.(signal)(halfCell1:cellLength));
    % Determines whether cell structure has to be flipped depending on signal
    % distribution and value given to the pole input.
               if (bottomFluo < upperFluo && pole == 0) || (bottomFluo > upperFluo && pole == 1)
                   cellToFlip = [cellToFlip cell]; %#ok<AGROW>
                  for frame = 1:lastFrame                            
    % Applies reorient function         
                       str = CL_getCell(cell,frame,cellList);
                       str2 = reorient(str);
                       cellList = CL_addCell(cell,frame,str2,cellList);
                   end
               end
             end
        end
    end
else
    for frames = 1:length(cellList.meshData)
        [~,cellIds] = CL_getFrame(frames,cellList);
        for cells = cellIds
            cellStructure = CL_getCell(cells,frames,cellList);
             % Calculates fluorescence intensity for input signal in each half of the cell. Only cells with
             % signal information will be included in the new structure.
            if ~isempty (cellStructure.(signal)) && isfield(cellStructure,'lengthvector')
                 cellLength = length(cellStructure.lengthvector);
                 halfCell = round(cellLength/2);
                 halfCell1 = halfCell+1;
                 bottomFluo = sum(cellStructure.(signal)(1:halfCell));
                 upperFluo = sum(cellStructure.(signal)(halfCell1:cellLength));
                  % Determines whether cell structure has to be flipped depending on signal
                  % distribution and value given to the pole input.
               if (bottomFluo < upperFluo && pole == 0) || (bottomFluo > upperFluo && pole == 1)
                     % Applies reorient function         
                     str2 = reorient(cellStructure);
                     cellList = CL_addCell(cells,frames,str2,cellList);
               end
            end
        end
    end
end
   
  
function str2 = reorient(str)
% this function reorients the data corresponding to one cell on a 
% single frame (without changing the polarity variable)
% 
names = fieldnames(str);
% flip the orientation of the data arrays (of every vertical array)
for i=1:length(names)
    t=['str2.' names{i} '=flipud(str.' names{i} ');'];
    eval(t);
end
% change the models for particular algorithms
if ismember('model',names) && ismember('algorithm',names)
   if ismember(str.algorithm,[1 2])
      str2.model(1) = mod(str.model(1) + pi,2*pi);
   elseif ismember(str.algorithm,3)
       str2.model(3) = mod(str.model(3) + pi,2*pi);
   elseif ismember(str.algorithm,4)
        s = size(str.model,1);
        s2 = floor(s/2);
        if s>1
           str2.model = [str.model(s2+1:s,:);str.model(1:s2,:)];
        end
    end
end
% reorient spots data (assuming the spots structure name starts with
% 'spots' and all spots sub-arrays are horizontal)
for i=1:length(names)
    if length(names{i})>=5 && strcmp('spots',names{i}(1:5))
       names2 = {};
       eval(['names2 = fieldnames(str.' names{i} ');'])
       if ismember('positions',names2)
          eval(['str2.' names{i} '.l = str.length-str.' names{i} '.l;'])
          eval(['str2.' names{i} '.positions = size(str.mesh,1)+1-str.' names{i} '.positions;'])
          eval(['names2 = fieldnames(str.' names{i} ');'])
          for j=1:length(names2)
              t=['str2.' names{i} '.' names2{j} '=fliplr(str2.' names{i} '.' names2{j} ');'];
              eval(t);
          end
        end
     end
end
% reorient lengthvector (the only aray for which flip is insufficient)
if ismember('lengthvector',names) && ismember('length',names)
   str2.lengthvector = str2.length-str2.lengthvector;
end
end

end