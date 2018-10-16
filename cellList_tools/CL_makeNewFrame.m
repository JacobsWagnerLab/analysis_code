function CL = CL_makeNewFrame(frame, CL)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function CL = CL_makeNewFrame(frame, CL)
%microbeTracker.v0.2.8
%@author:  Ahmad J Paintdakhi
%@date:    March 21, 2013
%@copyright 2012-2013 Yale University
%==========================================================================
%**********output********:
%CL:    A structure containing two fields meshData and cellId    
%**********Input********:
%CL:    A structure containing two fields meshData and cellId
%frame:  frame to be added to the cellList.
%==========================================================================
%This function adds an empty frame to the cellList if the frame is not
%in the cellList array.
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
lengthOfCellList = CL_getLength(CL);
%if frame already present skip the process and return;
if lengthOfCellList >= frame
   disp(['CL_makeNewFrame: ' num2str(frame) ' already exists!']);
   return;
end
%if frame to be added not in cellList, then add frame to cellList.
CL.meshData{frame} = {};
CL.cellId{frame} = [];
    
end