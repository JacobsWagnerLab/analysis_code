%% Plot the ensemble-averaged MSD.
%
% -About-
%   The method plotMeanMSD() plots the ensemble-averaged MSD in a loglog
%   scale. The slope of the loglog MSD plot indicates the mode of the
%   diffusion (e.g. normal, sub-diffusive, super-diffusive). The
%   y-intercept also scales with the population-averaged diffusion
%   coefficient. The resulted plot is a scatter plot and the color and
%   marker size of the scatter points can be specified by keywords in the
%   input. The default color is black, and the default marker size is 7.
% 
% -Input-
%   - obj: spt object
% 
% -Varargin-
%   - color: color used to plot the ensemble-averaged MSD
%   - markerSize: size for the markers
% 
% -Output-
%   A plot of the ensemble-averaged mean squared displacement in loglog
%   scales.
% 
% -Example-
%   % Plot the ensemble-averaged mean squared displacements
%   myParticle.plotMeanMSD();
%   % Plot with red scatter point with a marker size of 10
%   myParticle.plotMeanMSD('color','r','markerSize',10);
% 
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University


function plotMeanMSD(obj,varargin)

% Check if the ensemble-averaged MSD was calculated beforehand. If not
% calculate now.
if isempty(obj.meanMSD)
    obj.getMeanMSD();
end

% Default line color and marker size
color = 'k';
markersize = 7;

% Overwrite the default scatter styles based on the user input
for ii = 1: 2: length(varargin)
    switch lower(varargin{ii})
        case 'color'
            color = varargin{ii+1};
        case 'markersize'
            markersize = varargin{ii+1};
        otherwise
            error('Invalid argument. Valid options: color|markersize.');
    end
end

% Plot the ensemble-averaged MSD using physical units
plot(obj.meanMSD(:,1).*obj.frameTime,...
    obj.meanMSD(:,2).*(obj.pixelLength)^2,...
    's','color',color,'markersize',markersize); 

axis square;
grid on;
xlabel('Time [sec]');
ylabel('<MSD> [µm^2]');
set(gca,'xscale','log','yscale','log','linewidth',1.5,'fontsize',14);
end