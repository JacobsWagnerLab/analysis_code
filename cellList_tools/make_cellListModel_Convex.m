function cellList=make_cellListModel_Convex(cellList)
% cellList=make_cellListModel_Convex(cellList)
%   The function converts cell contours cpecified in the
%   cellList{..}{..}.model into convex polygons and saves a new model. 
% convinient if one wants convex cell contours instead of convoluted ones
% (for exampl created by pixel-based algorithms)
%==========================================================================
%
% INPUT:
% cellList - cellList (aka mesh) structure returned by Oufti
% OUTPUT:
% new_cellList - cellList with modified models.
%  NOTE: if model field is empty than nothing happens for that cell
%
%=========================================================================
%@author: Ivan Surovtsev  
%@date: 07.31.2015   
%@copyright 2012-2015 Yale University
%=========================================================================


n_frames=length(cellList.meshData);

for frame=1:n_frames

  MeshData=cellList.meshData{frame};
   IDs=cellList.cellId{frame};
  for cell=1:length(IDs)
    if ~isempty(MeshData{cell}.model)  
      ind=convhull(double(MeshData{cell}.model(:,1)),double(MeshData{cell}.model(:,2)));
      MeshData{cell}.model=[MeshData{cell}.model(ind,1),MeshData{cell}.model(ind,2)];
    end
  end

  cellList.meshData{frame}=MeshData;
  
end


end