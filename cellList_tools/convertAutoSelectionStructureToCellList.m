function cellList = convertAutoSelectionStructureToCellList(dataStructureFromAutoSelection,sizeOfImage)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function cellList = convertAutoSelectionStructureToCellList(dataStructureFromAutoSelection,sizeOfImage)
%@author:  Ahmad Paintdakhi
%@date:    October 01, 2014
%@copyright 2012-2015 Yale University
%==========================================================================
%**********output********:
%cellList:  data structure containing meshData and cellId fields used in
%Oufti application.
%**********Input********:
%dataStructureFromAutoSelection:  cell array structure gathered from 
%autoSelection function written by Manuel Campos.  Please see autoSelection
%function for more information.
%=========================================================================
% PURPOSE:  The purpose of this function is to convert a cell array of 
%structure back to cellList format.  The cell array of structures is formed
%after processing text formatted data of Oufti from high-throughput mode through
%auto selection procedure where a user selects only cells that go through a full 
%cell cycle after going through stringent conditions.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
count = 1; %counter for cells

%create a parse cellList array with fields meshData and cellId.  The length
%of the cellList is equal to the max number of frame found in input array.
cellList.meshData = cell(1,max(cat(2,dataStructureFromAutoSelection.frames)));
cellList.cellId = cell(1,max(cat(2,dataStructureFromAutoSelection.frames)));

%find ids of cells in the structure array
ids = cat(1,dataStructureFromAutoSelection.id);

%loop through all the frames and ids and start filling the cellList so that
%its format is usable in oufti.  For example, after changing the cell array 
%of structures to cellList format one will be able to use spot detection or 
%other modules in oufti on the data.
for frame = 1:length(cellList.meshData)
    for ii = 1:length(ids)
   
        dd = cat(1,dataStructureFromAutoSelection(ii).frames);
        xx = frame == dd;
        if sum(xx) == 0,  continue;end
        cellList.meshData{frame}{count}.mesh = dataStructureFromAutoSelection(ii).meshes{xx};
        cellList.meshData{frame}{count}.ancestors = dataStructureFromAutoSelection(ii).ancestor;
        cellList.meshData{frame}{count}.descendants = dataStructureFromAutoSelection(ii).progeny;
        cellList.meshData{frame}{count}.model = [cellList.meshData{frame}{count}.mesh(:,1:2);flipud(cellList.meshData{frame}{count}.mesh(2:end-1,3:4))];
        roiBox(1:2) = round(max(min(cellList.meshData{frame}{count}.model(:,1:2))-25,1));
        roiBox(3:4) = min(round(max(cellList.meshData{frame}{count}.model(:,1:2))+25),...
                     [sizeOfImage(2) sizeOfImage(1)])-roiBox(1:2);
        cellList.meshData{frame}{count}.box = roiBox;
        cellList.cellId{frame} = [cellList.cellId{frame} ids(ii)];
        cellList.meshData{frame}{count}.divisions = [];
        count = count + 1;
    end   
    count = 1;    
end

end