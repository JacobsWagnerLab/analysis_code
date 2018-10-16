function lengthFrame = CL_getFrameLength(frame, CL)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function lengthFrame = CL_getFrameLength(frame, CL)
%microbeTracker.v0.2.8
%@author:  Ahmad J Paintdakhi
%@date:    February 11, 2013
%@copyright 2012-2013 Yale University
%==========================================================================
%**********output********:
%lengthFrame:  length of a frame
%**********Input********:
%CL:        A structure containing two fields meshData and cellId
%frame:     frame number
%==========================================================================
%returns the length of a frame given by the input frame number
%-------------------------------------------------------------------------- 
%--------------------------------------------------------------------------  
lengthFrame = 0;
if ~CL_isEmptyFrame(frame, CL)
   lengthFrame = length(CL.meshData{frame});
end
end %function lengthFrame = CL_getFrameLength(frame, CL)
