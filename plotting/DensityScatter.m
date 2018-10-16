function figure_handle = DensityScatter(x,y,radius,x_name,y_name,method,contour_value,varargin)

%{
-About-
DensityScatter takes in two dimensions of data (x,y) and generates a
density-colored scatter plot according to the specifics defined by the
user. This function is an adaptation and combination of older versions of Manuel Campos's
denScatter functions.

-Inputs-
x:          input 1-d vector of some dimension
y:          y-values corresponding to x, a 1-d vector of the same length as x
radius:     positive number (does not have to be an integer) that defines the root with which the initial density
            estimates are transformed. Decrease for larger dark areas, Increase for
            larger light areas
x_name:     string specifying x-axis label
y_name:     string specifying y-axis label
method:     string specifying binning method to be added to the density
            scatter plot which may take on the following values
    'bin':      bins the data by grouping along x
    'nearest'	bins the data with a nearest neighbor approach. For each 
                point in x, N nearest points will be considered in the bin.
    'none':     simply plots the density scatter plot as is, without any
                additional items
contour:    integer defining type of isocontours displayed on graph.
            When selected, three isocontours defining regions of given
            densities will be displayed.
    0 : no isocontours
    1 : isocontours corresponding to 50-75-95 % probability of the data
    2 : isocontours corresponding to 1-2-3 standard deviations away from
    the highest density region

-varargin-
    'bins':       if scalar, the number of bins to group the data into. 
                (method == 'bin')

    'x_min':      minimum x-value used as boundary for binning.
                (method == 'bin')

    'x_max':      minimum x-value used as boundary for binning.
                (method == 'bin')

    'bin_method': Function handle used to calculate the value for each bin. If
                not provided will default to mean.
                available if (method == 'bin' || 'nearest')

    'error_method':   function handle used to calculate the error for each bin.
                    If not provided will default to standard deviation.
                    available if (method == 'bin')
 
    'n_neighbors':    The number of closest data points to use in binning 
                    available if (method == 'nearest')


-Outputs-
figure_handle:     handle to the generated figure.

-Example-

x = normrnd(2,1,5000,1);
y = normrnd(2,1,5000,1);

figure;
subplot(1,3,1)
DensityScatter(x,y,0.75,'X','Y','bin',0,'bins',10,'x_min',0,'x_max',4,'bin_method',@mean,'error_method',@std);
subplot(1,3,2)
DensityScatter(x,y,0.75,'X','Y','nearest',0,'n_neighbors',100,'bin_method',@median,'error_method',@mad);
subplot(1,3,3)
DensityScatter(x,y,0.75,'X','Y','none',1);

-Keywords-

Denisty scatter, plot, visualization

-Dependencies-

kde2dv2
BinDataInterval

-References-

Z. I. Botev, J. F. Grotowski, and D. P. Kroese. Kernel density estimation via diffusion. Ann Stat, 586
38(5):2916{2957, 2010.

-Author-
Sander Govers, 2017 November 6
%}

%% Establish default parameters

x_min=min(x);
x_max=max(x);
bins=10;
bin_method = @mean;
error_method = @std;
n_neighbors = 100;


%% Parse inputs

%choose the binning method to use based on the user's input
switch method
    case 'bin'
        visual_type = 'bin';
    case 'nearest'
        visual_type = 'nearest neighbors';
    case 'none'
        visual_type = 'none';        
end

%set parameters for selected binning method
for k = 1:length(varargin)
    
    if strcmpi(varargin{k}, 'bins')
        bins = varargin{k+1};
        if length(bins) ~= 1
            error('argument for bins was incorrectly formed')        
        end

    elseif strcmpi(varargin{k}, 'x_min')
        x_min = varargin{k+1};
        if length(x_min) ~= 1
            error('argument for x_min was incorrectly formed')        
        end

    elseif strcmpi(varargin{k}, 'x_max')
        x_max = varargin{k+1};
        if length(x_max) ~= 1
            error('argument for x_max was incorrectly formed')        
        end
        
    elseif strcmpi(varargin{k}, 'bin_method')
        if strcmp(class(varargin{k+1}),'function_handle')
            bin_method = varargin{k+1};
        else
            error('argument bin_method is not a callable function')
        end
    elseif strcmpi(varargin{k}, 'error_method')
        if strcmp(class(varargin{k+1}),'function_handle')
            error_method = varargin{k+1};
        else
            error('argument error_method is not a callable function')
        end
    elseif strcmpi(varargin{k}, 'n_neighbors')
        n_neighbors = varargin{k+1};
    end
end

%% Process data for generating scatter plot

%Clean up X and Y, getting rid of any NaN values
ixX=isnan(x);
ixY=isnan(y);
deleteRow=bsxfun(@or,ixX,ixY);
ix=deleteRow==0;

x=x(ix);
y=y(ix);

%Check wether input vectors are columns
if ~iscolumn(x)
    x=x';
end
if ~iscolumn(y)
    y=y';
end

%Calculation of data point density using kde2dv2 function (Manuel Campos,
%adapted from Zdravko Botev) in grid defined by a power of 2
[~,density,meshgrid_x,meshgrid_y]=kde2dv2([x y],2^8);

%Find respective coordinates of data points in meshgrid
[~,~,bin_idx_x] = histcounts(x,meshgrid_x(1,:));
[~,~,bin_idx_y] = histcounts(y,meshgrid_y(:,1));

%Use location of points to assign a color that corresponds to the density of
%the grid these data points occupy within the meshgrid
color_group=zeros(length(bin_idx_x),1);
for ii = 1:length(bin_idx_x)
    color_group(ii)=density(bin_idx_x(ii),bin_idx_y(ii));
end

%Reduce the number of color groups for increased speed during plotting
color_group=round(color_group*length(color_group)/3);

%Normalize distance
color_group=(color_group-min(color_group))/(max(color_group)-min(color_group));

%Adjust grouping according to radius
color_group=color_group.^(radius);
color_group=(color_group-min(color_group))/(max(color_group)-min(color_group));

%Define colormap
clrmap=flipud(([(0:0.009:0.9)',(0:0.009:0.9)',(0:0.009:0.9)']));

colors=clrmap(round(color_group*100)+1,:);


%% Generate scatter plot

figure_handle = scatter(x,y,15,colors,'filled');
set(gca,'fontname','arial','fontsize',16,'xcolor','k','ycolor','k');box on;
xlabel(x_name,'fontname','arial','fontsize',16,'fontweight','b','color','k');
ylabel(y_name,'fontname','arial','fontsize',16,'fontweight','b','color','k');
colormap(clrmap);
h=colorbar('YTick',[0 1],'fontname','arial','fontsize',14);
ylabel(h,'Relative density','fontsize',13,'fontweight','b','color','k');
hold on;


%% Determine bin type and plot corresponding errorbar
switch visual_type
    
    case 'bin'
        
        [x_out,y_out,error_out]=BinDataInterval(x,y,x_min,x_max,'bin','bins',bins,'bin_method',bin_method,'error_method',error_method);
        errorbar(x_out,y_out,error_out,'o','linewidth',1,'color',[1 0 0],'markerfacecolor',[1 0 0],'markersize',4);
        
    case 'nearest neighbors'
        
        [x_out,y_out,~]=BinDataInterval(x,y,x_min,x_max,'nearest','n_neighbors',n_neighbors,'bin_method',bin_method,'error_method',error_method);
        plot(x_out,y_out,'r','linewidth',1);
            
end

if contour_value~=0
    % Sort densities and calculate cumulative sum
    sorted_density = sort(density(:),'descend');
    cumulative_density=cumsum(sorted_density);
    
    if contour_value == 1
        %Find indexes corresponding to 50-75-95 % density probability of
        %the data
        ind1s = find(cumulative_density>=0.50*max(cumulative_density),1,'first');
        ind2s = find(cumulative_density>=0.75*max(cumulative_density),1,'first');
        ind3s = find(cumulative_density>=0.95*max(cumulative_density),1,'first');
        %Generate isocontours
        contour(meshgrid_x,meshgrid_y,density,[sorted_density(ind3s),sorted_density(ind2s), sorted_density(ind1s)],'--r');
    elseif contour_value == 2
        %Find indexes corresponding to 1-2-3 standard deviations (assuming normal distribution)
        ind1s = find(cumulative_density>=0.68268*max(cumulative_density),1,'first');
        ind2s = find(cumulative_density>=0.95449*max(cumulative_density),1,'first');
        ind3s = find(cumulative_density>=0.99730*max(cumulative_density),1,'first');
        contour(meshgrid_x,meshgrid_y,density,[sorted_density(ind3s),sorted_density(ind2s), sorted_density(ind1s)],'--r');
    else
        error('incorrect contour value')
    end
end

end