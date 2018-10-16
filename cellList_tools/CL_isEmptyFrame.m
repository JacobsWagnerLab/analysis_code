function isFrameEmpty = CL_isEmptyFrame(frame, CL)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function isFrameEmpty = CL_isEmptyFrame(frame, CL)
%microbeTracker.v0.2.8
%@author:  Ahmad J Paintdakhi
%@date:    March 22, 2013
%@copyright 2012-2013 Yale University
%==========================================================================
%**********output********:
%isFrameEmpty:  1 for true and 0 for false
%**********input********:
%CL:        A structure containing two fields meshData and cellId
%frame:     frame number
%==========================================================================
%The function checks if a frame is empty or not.
%-------------------------------------------------------------------------- 
%--------------------------------------------------------------------------
isFrameEmpty = 0;
if CL_isFrame(frame, CL)
    isFrameEmpty = isempty(CL.meshData{frame});
end

end