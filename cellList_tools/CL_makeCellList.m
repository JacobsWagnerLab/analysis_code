function newCellListFormat = CL_makeCellList(oldCellList)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%!!!!!!!!!!!!!!!!!!DO NOT USE THIS FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!DO NOT USE THIS FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!DO NOT USE THIS FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!DO NOT USE THIS FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!DO NOT USE THIS FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!DO NOT USE THIS FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!DO NOT USE THIS FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!DO NOT USE THIS FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!
%function newCellListFormat = CL_makeCellList(oldCellList)
%microbeTracker.v0.2.8
%@author:  Ahmad J Paintdakhi
%@date:    March 22, 2013
%@copyright 2012-2013 Yale University
%==========================================================================
%**********output********:
%newCellListFormat:  new cellList format containing fields meshData and
%cellId.
%**********input********:
%oldCellList:	old cellList as a cell array.
%==========================================================================
%The function makes the new cellList format from the old one with fields
%meshData and cellId.
%-------------------------------------------------------------------------- 
%--------------------------------------------------------------------------

%if oldCellList is empty then make new cellList with empty cell fields and
%return without invoking further lines in the function.
if isempty(oldCellList),newCellListFormat.meshData = {[]}; newCellListFormat.cellId = {[]}; return; end

newCellListFormat.meshData = cell(1,length(oldCellList));
newCellListFormat.cellId   = cell(1,length(oldCellList));
for frame = 1:length(nonzeros(~cellfun(@isempty,oldCellList))')
    I = ~cellfun('isempty', oldCellList{frame}); % Non-empty elements.
    c    = oldCellList{frame}(I);
    nr   = 1:length(I);
    ids  = nr(I);
    newCellListFormat = CL_addFrame(frame, c, ids, newCellListFormat);
end
end