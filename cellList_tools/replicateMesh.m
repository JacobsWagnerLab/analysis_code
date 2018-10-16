function new_cellList=replicateMesh(cellList,n_frames,option)
%
%--------------------------------------------------------------------------
% function []=replicateMesh(n_frames)
%
% Replicate meshes provided in cellList for a given number of frames n_frames. 
% for example for FRAP-analysis, when phase image is taken only before FRAP, 
% or before and after FRAP, and one need to generate a cellList covering FRAP experiment
% Assumes "oufti" cellList fromat. If you have older ('microbeTracker')
% format, convert it by using CL_makeCellListNew function from Oufti
% package.
%==========================================================================
%
% INPUT:
% cellList - cellList (aka mesh) structure returned by Oufti
% n_frames - number of frames in new cellList
% option   - optional, if option=1 then only first frame of celllist will be used for replication 
%
% OUTPUT:
% new_cellList - replicated cellList:
%   if input cellList has 1 frame, this frame will be replicated n_frames times
%   if input cellList has 2 frames, they will be first and last frames,
%   if input has more then 2 frames only first 2 will be used, other will be ignored 
%   frames in between will have everyting as in frames except, .mesh and .model fileds
%   will be linearly interpolated between meshes in first and last frames
%
%=========================================================================
%@author: Ivan Surovtsev  
%@date: 04.09.2015   
%@copyright 2012-2015 Yale University
%=========================================================================
% Publication: 
% Publication Author(s):
% Publication Title:
%--------------------------------------------------------------------------

meshData=cellList.meshData;
cellId=cellList.cellId;

if nargin<3 % to switch to "only first frame" mpode if option is provided 
  n_frames_0=length(cellId);
else
  n_frames_0=option;  
end

for ii=1:n_frames
  new_meshData{ii}=meshData{1};
  new_cellId{ii}=cellId{1};
end
new_cellList.meshData=new_meshData;
 new_cellList.cellId=new_cellId;
 
switch n_frames_0
  
  case 1
    disp(['First frame has been replicated into ',num2str(n_frames),' frames in new cellLIst'])
  
  otherwise
        
    new_cellId_0=cellId{1}; new_cellId_1=cellId{2};
    if length(new_cellId_0)==length(new_cellId_1)
      check=(new_cellId_0==new_cellId_1);
      if sum(check)==length(check)
        for cell=1:length(new_cellId_0)
          model_0=meshData{1}{cell}.model; model_1=meshData{2}{cell}.model;  
          mesh_0=meshData{1}{cell}.mesh; mesh_1=meshData{2}{cell}.mesh;
          for frame=1:length(new_meshData)  
            new_meshData{frame}{cell}.model=model_0+(model_1-model_0)*(frame-1)/(n_frames-1);
            new_meshData{frame}{cell}.mesh=mesh_0+(mesh_1-mesh_0)*(frame-1)/(n_frames-1);
          end
        end
        if n_frames_0>2, disp(['New cellList of ', num2str(n_frames),' frames have been created by interpolation between first and second frames of input cellList']); end   
      else
        disp('Not the same cells appears at first and last frames - simply replicated first frame')  
        return  
      end
    else
      disp('Not the same cells appears at first and last frames - simply replicated first frame')  
      return 
    end
end
    
        

end