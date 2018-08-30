%% Plot individual MSD curves.
%
% -About-
%   The method plotMSD() plots the individual MSD curves in a loglog scale.
%   The user can specify the indices of the MSDs to be plotted and the plot
%   styles using keywords. The default plot style results a jet-colormap
%   colored dots without any line connecting.
% 
% -Input-
%   - obj: spt object
% 
% -Varargin-
%   - index: an index array used to specify which MSD curves to be plotted
%   - color: color used for all the MSD curves
%   - marker: marker style for the MSD curves
%   - markerSize: size for the markers used above
% 
% -Output-
%   A plot of individual MSD curves in loglog scales.
%
% -Example-
%   % Plot all the MSD curves
%   myParticle.plotMSD();
%   % Plot the first 30 MSD curves
%   myParticle.plotMSD('index',1:30);
%   % Plot the first 10 MSD curves using red square markers
%   myParticle.plotMSD('index',1:10,'color','r','marker','s');
% 
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function plotMSD(obj,varargin)

% Check if the MSD are calculated beforehand. If not, calculate now.
if isempty(obj.MSD)
    obj.getMSD();
end

% Default parameters
% Indicies of the MSDs to be plotted
index = 1:size(obj.MSD,1);
% Plot styles
% If multiple tracks plotted, use the jet colormap. Else, use blue.
% marker = '.';
% markersize = 7;
% Memory for color to be determined later
color = '';

% Overwrite the default parameters based on the user input
for ii = 1 :2: length(varargin)
    switch lower(varargin{ii})
        case 'color'
            color = varargin{ii+1};
%         case 'markersize'
%             markersize = varargin{ii+1};
%         case 'marker'
%             marker = varargin{ii+1};
        case 'index'
            index = varargin{ii+1};
        otherwise
            error('Invalid argument. Valid options: index|color|marker|markersize.');
    end
end

% Depending on the final color set, use either a jet colormap or a single
% color
if (isempty(color)) && (length(index) == 1)
    color = repmat([0 0 0],[length(index),1]);
elseif (isempty(color)) && (length(index) ~= 1)
    color = jet(length(index));
else
    switch color
        case 'k'
            color = [0 0 0];
        case 'y'
            color = [1 1 0];
        case 'm'
            color = [1 0 1];
        case 'c'
            color = [0 1 1];
        case 'r'
            color = [1 0 0];
        case 'g'
            color = [0 1 0];
        case 'b'
            color = [0 0 1];
        otherwise
            color = [0 0 0];
    end
    color = repmat(color,[length(index),1]);
end
    

% Go through the indexed MSD and plot each of them using the specified plot
% styles.

for ii = 1:length(index)
    % Convert to physical units
    delay = obj.MSD{index(ii)}(:,1) .* obj.frameTime;
    msd = obj.MSD{index(ii)}(:,2) .* (obj.pixelLength)^2;
 
    % Plot the MSD curves
    plot(delay,msd,'color',[color(ii,:),0.1]);
    %plot(delay,msd,'color',color(ii,:),'linestyle','none','marker',marker,'markersize',markersize);
    hold on;
end


axis square;
grid on;
xlabel('Time [sec]');
ylabel('MSD [µm^2]');
set(gca,'xscale','log','yscale','log','linewidth',1.5,'fontsize',14);
end
