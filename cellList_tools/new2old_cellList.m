function cellList_old=new2old_cellList(cellList_new)

% this function converts current cellList structure (aka "new cellList") into
% old cellList structure (aka "Oleksii's cellList"). 
% Specifically:
%   - cellList_old has only data by frames, which are data by cells 
%   - converts single-type variable into double-variable
%   - add extra filed, such as ".lengthvector" etc
%   - add empty entries so cell numbering match cell ID
%
% Version: 03-12-2013
%
% by Ivan Surovtsev

cellList=cellList_new.meshData;
cell_ID=cellList_new.cellId;

CL_fields={'model','mesh','box','signal0','signal1','signal2'};

for ii=1:length(cellList)
    
    for jj=1:length(cellList{ii})
        kk=cell_ID{ii}(jj);
        cellList_old{ii}{kk}=cellList{ii}{jj};
        for aa=CL_fields; 
            if isfield(cellList{ii}{jj},aa{1})
                 cellList_old{ii}{kk}.(aa{1})=double(cellList{ii}{jj}.(aa{1}));
            end
        end
        cellList_old{ii}{kk}=getextradata(cellList_old{ii}{kk});
        
    end    
end
end

    
