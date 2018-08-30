function [xy_base_error_poly_basegons, X, y_base] = Error_Fill_Pr_Hist(indat,bins)
%{
-About-
Calculate a histogram with area of 1 with Poisson counting error.

-Inputs-
indat:  input data for the histogram calculation
bins:  the number of histogram bins to be produced. If bins is not
specified, the number of bins will be sqrt(length of data)

-varargin-
N/A

-Outputs-
xy_base_error_poly:     a cell array of polygons to be plotted. Each
polygon spans the high error and low error for each bin and is meant to be
plotted either in a loop or with cellfun, see the below example.
X:  The histogram X data
y_base: the histogram y data

-Example-
dx = randn(1,100000);
[err_poly, ~, ~] = eFillStepPr((dx),120);
%since there will likely be multiple polygons, use the 'hold on' command
figure,hold on
color = [0.2,0.2,0.7];
cellfun(@(x) fill(x(:,1),x(:,2),color,'EdgeColor',color), err_poly)
xlabel('Values')
ylabel('Frequency')
set(gca,'yscale','log','fontsize',13,'xtick',[-5:5])
axis tight
   
-Supplementary-
N/A

-Keywords-
histogram, Poisson error

-Dependencies-
None

-References-
Used to generate histogram plots in figures 2-5 of Parry et al, Cell, 2014

-Author-
Brad Parry
%}


if nargin == 1
    bins = floor(sqrt(length(indat)));
end

%compute base histogram data
[y_base, X] = hist(indat, bins);

%prepare a few anony_basemous helper functions for the integrattion and error
%calculations
Area_Int_One = @(x,y_base,hdata) y_base/sum(y_base.*diff([min(hdata) x(1:end - 1) + diff(x)/2 max(hdata)]));
Area_Int_One_H = @(x,y_base,hdata) (y_base + sqrt(y_base))/sum(y_base.*diff([min(hdata) x(1:end - 1) + diff(x)/2 max(hdata)]));
Area_Int_One_L = @(x,y_base,hdata) (y_base - sqrt(y_base))/sum(y_base.*diff([min(hdata) x(1:end - 1) + diff(x)/2 max(hdata)]));

y_base_high = Area_Int_One_H(X,y_base,indat);
y_base_low = Area_Int_One_L(X,y_base,indat);
y_base = Area_Int_One(X,y_base,indat);

contigs = [0 (y_base ~= 0) 0];
start_inds = strfind(contigs, [0, 1]);
end_inds = strfind(contigs, [1, 0]) - 1;

xy_base_error_poly_basegons = cell(1,length(start_inds));
for q = 1:length(start_inds)
    roi_y_base = [y_base_high(start_inds(q):end_inds(q)), y_base_low(end_inds(q):-1:start_inds(q))]';
    roi_x = [X(start_inds(q):end_inds(q)), X(end_inds(q):-1:start_inds(q))]';
    roi_y_base = cat(1, roi_y_base, roi_y_base(1)); 
    roi_x = cat(1, roi_x, roi_x(1));
    roi_x(roi_y_base == 0) = [];
    roi_y_base(roi_y_base == 0) = [];
    xy_base_error_poly_basegons{q}(1:size(roi_x,1),1:2) = [roi_x'; roi_y_base']'; 
end