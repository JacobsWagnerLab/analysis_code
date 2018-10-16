function slicedCellList = CL_sliceFrames(frameRange, CL)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function slicedCellList = CL_sliceFrames(frameRange, CL)
%microbeTracker.v0.2.8
%@author:  Ahmad J Paintdakhi
%@date:    March 27, 2013
%@copyright 2012-2013 Yale University
%==========================================================================
%**********output********:
%slicedCellList:    A section of a cellList containing two fields meshData and cellId    
%**********Input********:
%CL:    A structure containing two fields meshData and cellId
%frameRange:  frame range to be used for sectioning such as [1 10], which
%means data from frame 1-10 is desired.
%==========================================================================
%The function returns a section of cellList given by frame range.
%-------------------------------------------------------------------------- 
%--------------------------------------------------------------------------

if CL_isFrame(frameRange(1), CL) && CL_isFrame(frameRange(2), CL);
   slicedCellList.meshData   = CL.meshData(frameRange(1):frameRange(2));
   slicedCellList.cellId = CL.cellId(frameRange(1):frameRange(2));
else
   slicedCellList.meshData    = {};
   slicedCellList.cellId = {};
end

end

