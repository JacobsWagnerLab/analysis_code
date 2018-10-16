function [div_num,f_div1]=is_any_div2(cellList,cell,f_in)
% to check if there any cell division occured for a specified cell 
%
% INPUT
% cellList = cellList structure
% cell = cell number, aka cell ID, to check
% f_in = starting frame, runs from this rame to the final frame of
%        cellList, if f_in <1 then starts from frame=1
% OUTPUT
% div_num = number of divisions since f_in to the last frame, 0 if no divisions)
% f_div1 = frame at which first division after frame=f_in  occured
%
% Version: 12-20-11
%
% by Ivan Surovtsev
  
  div_num=0;
  f_fin=length(cellList);
  f_div1=f_fin+1;
  
  [out,fr_1]=is_cell2(cellList,cell);
  
  if out
      
    f_in=max([1,f_in,fr_1]);
    
    if ~isempty(cellList{f_in})
      div_0=length(cellList{f_in}{cell}.divisions);
    else
      while isempty(cellList{f_in})
          f_in=f_in+1;
      end
      div_0=length(cellList{f_in}{cell}.divisions);
    end
   div_t=div_0;
    div_num=div_t-div_0;
      
    for frame=f_in:f_fin    
        if ~(cell>length(cellList{frame}))
          if ~isempty(cellList{frame}{cell}) && length(cellList{frame}{cell}.mesh)>1 
            div=length(cellList{frame}{cell}.divisions);
            if div>div_t
              div_t=div;
              f_div1=cellList{frame}{cell}.divisions(div_0+1);
              div_num=div_t-div_0;   
            end
          end
        end
    end
     
    if div_num==0
      f_div1=f_fin+1;
    end
    
  else
      
    return
  
  end
      

  %end
  
end