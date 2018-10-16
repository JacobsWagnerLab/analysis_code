function [spc_max,f_max]=is_any_spots2(cellList,cell,f_in)
% to check if there are any spots AND meshes for a specified cell
%
% INPUT
% cellList = cellList structure
% cell = cell number, aka cell ID, to check
% f_in = starting frame, runs from this rame to the final frame of
%        cellList, if f_in <1 then starts from frame=1
% OUTPUT
% spc_max = maximum number of spots per cell within analyzed frames, 0 if
%          no spots
% f_max = frame at which first spc_max number of spots observed, 0 if no spots 
%
% Version: 12-20-11
%
% by Ivan Surovtsev

  
  spc_max=0;
   f_max=0;
  
  f_in=max(1,f_in);
   f_fin=length(cellList);
     
  for frame=f_in:f_fin    
      if ~(cell>length(cellList{frame}))
        if ~isempty(cellList{frame}{cell}) && length(cellList{frame}{cell}.mesh)>1 
          spc=length(cellList{frame}{cell}.spots.l);
          if spc>spc_max
            spc_max=spc; f_max=frame;  
          end
        end
      end
  end
  
  
end