function segintarray = CL_segHist(cellList)
%------------------------------------------------------------------
%------------------------------------------------------------------
% author: Sangjin Kim
% date: July 3, 2013
% copyright 2012-2013 Yale University
%======================================================
%***************output******************************
% segintarray = an array of all segment intensities
% figure, histogram of segment intensities
% number, max segment intensity, used later for spotThreshold
%***************input******************************
% cellList
%======================================================
%This script will give distribution of segment intensities in all cells in
%the cellList.
%
segintarray = [];
for frame=1:length(cellList)
    for cell=1:length(cellList{frame})
        if ~isempty(cellList{frame}{cell}) && length(cellList{frame}{cell}.mesh)>4
            segintarray = [segintarray reshape(cellList{frame}{cell}.signal2./cellList{frame}{cell}.steparea,1,[])];
        end
    end
end

c = 0:(max(segintarray)/sqrt(length(segintarray))):max(segintarray);
h1 = hist(segintarray,c);
figure, bar(c,h1')
xlabel('Mean intensity in a segment')
ylabel('%');

disp('max seg intensity');
disp(max(segintarray));
