function isContourAvailable = CL_isContourCell(cellId, frame, CL)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function isContourAvailable = CL_isContourCell(cellId, frame, CL)
%microbeTracker.v0.2.8
%@author:  Ahmad J Paintdakhi
%@date:    March 22, 2013
%@copyright 2012-2013 Yale University
%==========================================================================
%**********output********:
%isContourAvailable:  1 for true and 0 for false
%**********input********:
%CL:        A structure containing two fields meshData and cellId
%cellId:    id of a cell to be accessed.  The id is located in CL.cellId
%frame:     frame number
%==========================================================================
%The function first checks if a cell is available, if that is the case then
%it checks for the availability of a contour.
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
cellStructure =  CL_getCell(cellId, frame, CL);    
isContourAvailable = ~isempty(cellStructure) && (isfield(cellStructure,'contour') && length(cellStructure.contour)>1);
end