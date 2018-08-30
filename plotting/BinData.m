function [Xout, Yout] = BinData(x,y, method, varargin)
%{
-About-
BinData takes in two dimensions of data (x,y) and bins the data along x 
according to the method specificed by 'method' and parameteres passed. 
BinData returns the binned coordinates Xout and Yout

-Inputs-
x:          x-coordinates of input data, must be the same size as y
y:          y-coordinates of input data, must be the same size as x
method:     a string that specifies the type of binning to be performed. Method
            may take on the following values:
	'bin'		Bins the data by grouping along x in bins of a fixed size
                specified by the bin argument (see below)
	'nearest'	Bins the data with a nearest neighbor approach. For each 
                point in x, N nearest points will be considered in the bin.
                N is determined by the n_neighbors argument, see below.
	'gaussian'	Weighted average of data points where the weights are 
                computed as the distance along a Gaussian profile. The
                shape of the Gaussian distribution is controlled by sigma,
                see below.

-Varargin-
bins:       If scalar, the number of bins to group the data into. 
            If matrix, should be a 2 by n matrix where the first row 
            deliniates the left bounds of each bin and the second row 
            deliniates the right bound of each bin (method = 'bin' only)

bin_method: Function handle used to calculate the value for each bin. If
            not provided will default to mean.
            available if (method == 'bin' || 'nearest')

n_neighbors: The number of closest data points to use in binning 
            available if (method == 'nearest')

sigma:	  Controls the shape of the Gaussian function used to determine 
            data distances. available if (method == 'gaussian')

-Outputs-
Xout    the input data, x, binned according to the user-specified method
Yout    the input data, y, binned according to the user-specified method

-Example-
#to run the example, either call BinData() without any inputs, or copy and
#paste the following code in the matlab prompt

-Keywords-
plot, bin, group, nearest neighbors, similarity, gaussian binning

x = 1:.01:3*pi;
y = sin(x) + .7*randn(size(x));

figure
h = get(gcf,'position');
set(gcf,'position',[int32(h(1)/3), h(2), h(3)*2.5, h(4)])
%Gaussian binning with std = .1
[Xout, Yout] = BinData(x, y, 'gaussian','sigma',.2);
subplot(1,3,1)
plot(x,y,'.k')
hold on
plot(Xout,Yout,'m','linewidth',3)
title('Gaussian, sigma = 0.2')

%nearest neighbors
[Xout, Yout] = BinData(x, y, 'nearest','n_neighbors',100,'bin_method',@median);
subplot(1,3,2)
plot(x,y,'.k')
hold on
plot(Xout,Yout,'m','linewidth',3)
title('nearest, N = 100')

%bin in 50 bins
[Xout, Yout] = BinData(x, y, 'bin','bins',50,'bin_method',@mean);
subplot(1,3,3)
plot(x,y,'.k')
hold on
plot(Xout,Yout,'m','linewidth',3)  
title('bins, N = 50')      
   
-Author-
Brad Parry, Feb 11, 2016, modified June 22, 2017
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       establish default parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%method = bin
if nargin == 0
    run_example
    return
end
bins = 10;
bin_method = @mean;

%method = nearest neighbors
n_neighbors = round( length(x) / 10 );

%method = gaussian weight
sigma = .3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       parse inputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%choose the binning method to use based on the user's input
switch method
    case 'bin'
        binning_type = 'bin';
    case 'nearest'
        binning_type = 'nearest neighbors';
    case 'gaussian'
        binning_type = 'gaussian weight';
end

%set parameters for selected binning method
for k = 1:length(varargin)
    
    if strcmpi(varargin{k}, 'bins')
        bins = varargin{k+1};
        if length(bins) == 1
            [bins_L, bins_R] = ConstructBins(x, bins);
            
        elseif size(bins,1) == 2 && size(bins,2) > 1
            bins_L = bins(1,:);
            bins_R = bins(2,:);
            
        else
            error('argument for bins was incorrectly formed')
            
        end
        
    elseif strcmpi(varargin{k}, 'bin_method')
        if strcmp(class(varargin{k+1}),'function_handle')
            bin_method = varargin{k+1};
        else
            error('argument bin_method is not a callable function')
        end
        
    elseif strcmpi(varargin{k}, 'n_neighbors')
        n_neighbors = varargin{k+1};
        
    elseif strcmpi(varargin{k}, 'sigma')
        sigma = varargin{k+1};

    end
end

%all methods require sorted data: sort data by x dimension
[~,ix] = sort(x);
x = x(ix);
y = y(ix);

%Call the appropriate binning method
switch binning_type
    
    case 'bin'
        
        [Xout, Yout] = AggregateData(x, y, bins_L, bins_R, bin_method);
        
    case 'nearest neighbors'
        
        [Xout, Yout] = NearestNeighbors(x, y, n_neighbors, bin_method);
        
    case 'gaussian weight'
        
        [Xout, Yout] = GaussianWeightX(x, y, sigma);        
        
end

end

function [bins_L, bins_R] = ConstructBins(x,bins)
%construct the number of bins specified by the user along the data's x-axis
%the bins returned will be uniformly spaced along the data
bins = round(1 : (length(x)-1)/bins : length(x));
bins_L = bins(1:length(bins)-1);
bins_R = bins(2:end);
bins_R(1:end-1) = bins_R(1:end-1) - 1;
end

function [wx, wy] = AggregateData(x, y, bins_L, bins_R, bin_method)
% perform data aggregation specified by the anonymous function bin_method
% over each of the k bins in bins_L and bins_r. The k'th bin is specifed by
% the data lying in [bins_L(k), bins_r(k)]

%initialize
wx = zeros(1,length(bins_L));
wy = zeros(1,length(bins_R));
for k = 1:length(bins_L)
    
    %call the anonymous function passed to transform the data
    wx(k) = bin_method(x(bins_L(k):bins_R(k)));
    wy(k) = bin_method(y(bins_L(k):bins_R(k)));
    
end

end

function [wx, wy] = GaussianWeightX(x, y, sigma)
% GaussianWeightX constructs a weighted average of all points (x,y) with
% weights determined by the normal equation. The width of the weights are 
% specified by sigma. If samples in x are uniformly spaced, this routine 
% could be optimized by convolution with a Gaussian kernel. The method 
% below, however, is universal in that the spacing of the data does not 
% matter.

%prepare an anonymous function for convenience
tile = @(x) repmat(x(:), [1, length(x)]);

%make sure the data is the correct shape (i.e., column vectors)
x = x(:);
y = y(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% the next segment measures the point to point difference for all points in
% the data and then scales each batch of point-to-point differences by
% their standard deviation. The weight of all points for the calculation at
% a single point is then determined by N(0, sigma), where N returns the
% weights from the unnormalized Normal equation:
%
% P(x) = e^{{{ - \left( {scaled\_L1\_dist } \right)^2 } \mathord{\left/ 
% {\vphantom {{ - \left( {x - \mu } \right)^2 } {2\sigma ^2 }}} \right. 
% \kern-\nulldelimiterspace} {2\sigma ^2 }}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% measure the L1 distance from each point to all other points in the
% dataset
L1_dist = repmat(x(:), [1, length(x)]) - repmat(x(:)', [length(x), 1]);
% scale each column of distance by its standard deviation so that the
% point-to-point sigma can be measured correctly
scaled_dist = L1_dist ./ repmat(std(L1_dist, 1), [size(L1_dist,1), 1]);

% apply transformed distaneces to weight each point in a sigma neighborhood
% based on a 1-D Gaussian distribution
weights = exp( -scaled_dist.^2 ./ (2*sigma^2) );

% Use the Gaussian weights, weights, to calculate weighted averages
wx = sum(weights.*tile(x)) ./ sum(weights);
wy = sum(weights.*tile(y)) ./ sum(weights);
end

function [wx, wy] = NearestNeighbors(x, y, n_neighbors, bin_method)
% identify n_neighbors nearest neighbors for each data point and aggregate
% them with the method specified by the anonymous function bin_method

wx = zeros(size(x));
wy = zeros(size(y));

for k = 1:length(x)

    dx = (x(k) - x).^2;
    
    [~,ix] = sort(dx);
    ix = ix(1:n_neighbors);

    wx(k) = bin_method(x(ix));
    wy(k) = bin_method(y(ix));

end

end

function run_example()
x = 1:.01:3*pi;
y = sin(x) + .7*randn(size(x));

figure
h = get(gcf,'position');
set(gcf,'position',[int32(h(1)/3), h(2), h(3)*2.5, h(4)])
%Gaussian binning with std = .1
[Xout, Yout] = BinData(x, y, 'gaussian','sigma',.2);
subplot(1,3,1)
plot(x,y,'.k')
hold on
plot(Xout,Yout,'m','linewidth',3)
title('Gaussian, sigma = 0.2')

%nearest neighbors
[Xout, Yout] = BinData(x, y, 'nearest','n_neighbors',100,'bin_method',@median);
subplot(1,3,2)
plot(x,y,'.k')
hold on
plot(Xout,Yout,'m','linewidth',3)
title('nearest, N = 100')

%bin in 50 bins
[Xout, Yout] = BinData(x, y, 'bin','bins',50,'bin_method',@mean);
subplot(1,3,3)
plot(x,y,'.k')
hold on
plot(Xout,Yout,'m','linewidth',3)  
title('bins, N = 50')
end


