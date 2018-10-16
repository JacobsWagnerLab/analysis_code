function isNonEmptyFrame = CL_isNonEmptyFrame(frame, CL)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function isNonEmptyFrame = CL_isNonEmptyFrame(frame, CL)
%microbeTracker.v0.2.8
%@author:  Ahmad J Paintdakhi
%@date:    March 22, 2013
%@copyright 2012-2013 Yale University
%==========================================================================
%**********output********:
%isNonEmptyFrame:  1 for true and 0 for false
%**********input********:
%CL:        A structure containing two fields meshData and cellId
%frame:     frame number
%==========================================================================
%The function finds if the frame is empty or not.
%-------------------------------------------------------------------------- 
%--------------------------------------------------------------------------
isNonEmptyFrame = CL_isFrame(frame, CL) && ~CL_isEmptyFrame(frame, CL);
end