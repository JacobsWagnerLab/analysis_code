function [spc_max,f_max,frames]=is_any_spots3(cellList,cell,f_in)
% to check if there are any spots AND meshes for a specified cell
% returns arrays of frames where cell exists, where 1 spot detected, where 2 spots detected etc 
% frames.all= all frames in wich cell mesh exists
% frames.0= frames w 0 spots
% frames.1= frames w 1 spots
% frames.2= frames w 2 spots
% frames.3= frames w >2 spots
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
% Initial Version: 12-20-11
% Updated versions: 05.28.13 fixed bug: by error all frames with >=2 spots were count as 2-spots
%
% by Ivan Surovtsev

  
  spc_max=0;
   f_max=0;
  frames.all=[]; frames.spc0=[]; frames.spc1=[]; frames.spc2=[]; frames.spc3=[];
   
  f_in=max(1,f_in);
   f_fin=length(cellList);
     
  for frame=f_in:f_fin    
      if ~(cell>length(cellList{frame}))
        if ~isempty(cellList{frame}{cell}) && length(cellList{frame}{cell}.mesh)>1 && isfield(cellList{frame}{cell},'spots') 
          frames.all=[frames.all,frame];  
          spc=length(cellList{frame}{cell}.spots.l);
          if spc>spc_max
            spc_max=spc; f_max=frame;  
          end
          switch spc
              case 0
                  frames.spc0=[frames.spc0,frame];
              case 1
                  frames.spc1=[frames.spc1,frame];
              case 2
                  frames.spc2=[frames.spc2,frame];
              otherwise
                  frames.spc3=[frames.spc3,frame];
          end
        end
      end
  end
  
  
end