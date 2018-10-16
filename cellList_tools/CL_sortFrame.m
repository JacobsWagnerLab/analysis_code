function CL = CL_sortFrame(frame, CL)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function CL = CL_sortFrame(frame, CL)
%microbeTracker.v0.2.8
%@author:  Ahmad J Paintdakhi
%@date:    March 27, 2013
%@copyright 2012-2013 Yale University
%==========================================================================
%**********output********:
%CL:    A sorted frame is added to the already CL - cellList.    
%**********Input********:
%CL:    A structure containing two fields meshData and cellId
%==========================================================================
%The function sorts a cellList - CL only at a given frame number.
%-------------------------------------------------------------------------- 
%--------------------------------------------------------------------------

[meshData cellIds] = CL_getFrame(frame, CL);
[cellIds I] = sort(cellIds);
meshData = meshData(I);
CL = CL_addFrame(frame, meshData, cellIds, CL); % overwrites the data in old frame.
end