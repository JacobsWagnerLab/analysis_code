function [] = CL_spotThreshold(cellList,ch,threshold)
%------------------------------------------------------------------
%------------------------------------------------------------------
% author: Sangjin Kim
% date: July 3, 2013
% copyright 2012-2013 Yale University
%======================================================
%***************output******************************
% cellList1 = same as intput cellList, except some low intensity spots will
% be removed
%***************input******************************
% cellList (containing spot information)
% ch = 'signal1' or 'signal2'
% threshold = value for cutting spot intensity (obtained from CL_segHist)
% Change threshold value based on the max segment value
%======================================================
%This script will remove spots whose intensities are lower than threshold value.
% Change threshold value based on the max segment value

for frame=1:length(cellList)
    for cell=1:length(cellList{frame})
        if ~isempty(cellList{frame}{cell}) && isfield(cellList{frame}{cell},'spots') &&isfield(cellList{frame}{cell},ch) ...
                && ~isempty(cellList{frame}{cell}.spots.x)
            prf = cellList{frame}{cell}.signal1./cellList{frame}{cell}.steparea;
            for i=1:3
                prf = prf*0.5+prf([1 1:end-1])*0.25+prf([2:end end])*0.25; %smoothing
            end
            
            pos = cellList{frame}{cell}.spots.positions;
            prf = [0;prf];
            
            ind = prf(pos+1)>threshold; %1 or 0
             
            cellList{frame}{cell}.spots.positions = cellList{frame}{cell}.spots.positions(ind); %takes only 1s
            cellList{frame}{cell}.spots.x = cellList{frame}{cell}.spots.x(ind);
            cellList{frame}{cell}.spots.y = cellList{frame}{cell}.spots.y(ind);
            cellList{frame}{cell}.spots.l = cellList{frame}{cell}.spots.l(ind);
            cellList{frame}{cell}.spots.d = cellList{frame}{cell}.spots.d(ind);
            cellList{frame}{cell}.spots.magnitude = cellList{frame}{cell}.spots.magnitude(ind);
            cellList{frame}{cell}.spots.h = cellList{frame}{cell}.spots.h(ind);
            cellList{frame}{cell}.spots.w = cellList{frame}{cell}.spots.w(ind);
            cellList{frame}{cell}.spots.b = cellList{frame}{cell}.spots.b(ind);
        end
    end
end

cellList1 = cellList;