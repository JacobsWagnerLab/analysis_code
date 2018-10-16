function [proccells CL] = CL_packFirstFrame(proccells, CL)
    % Changes the cell ids to [1, 2, ...] and renumbers the list of cells in
    % proccells to follow the new numbering.

    [cells, ids] = CL_getFrame(1, CL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ahmad.P May 17 2012 keep old numbers.
%     newIds = 1:length(ids);
%     
%     nrCells = length(proccells);
%     %if nrCells > 1
%         for i=1:length(proccells)
%             proccells(i) = newIds(ids==proccells(i));
%         end
%     %else
%     %        proccells(i) = newIds(ids==proccells);
%     %end
%CL = CL_addFrame(1, cells, newIds, CL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    CL = CL_addFrame(1, cells, ids, CL);
end