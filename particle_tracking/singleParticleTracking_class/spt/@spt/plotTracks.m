%% Plot individual tracks.
%
% -About-
%   The method plotTracks() plots the particle trajectories stored in the
%   object. The user can provide an index array to specify which tracks to
%   be plotted.
% 
% -Input-
%   - obj: spt object
% 
% -Varargin-
%   An index array specifying which tracks to be plotted. If the index was
%   not provided, all tracks stored in the object are plotted.
% 
% -Output-
%   A plot of particle trajectories
% 
% -Example-
%   % Plot all tracks
%   myParticle.plotTracks();
%   % Plot the 3rd track
%   myParticle.plotTracks(3);
%   % Plot the first 10 tracks
%   myParticle.plotTracks(1:10);
% 
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University


function plotTracks(obj,varargin)

% Check if an index array for which tracks to be plotted are given
% If no index is given, plot all tracks.
if isempty(varargin)
    index = 1:length(obj.tracks);
else
    index = varargin{1};
end

% Color map for the tracks
% If multiple tracks plotted, use the jet colormap. Else, use blue.
if (length(index) == 1)
    color = 'b';
else
    color = jet(length(index));
end

% Plot all tracks
clf;
for ii = 1: length(index)
    
    % Extract a single track
    thisTrack = obj.tracks{index(ii)};
    
    % Convert to physical units
    x = thisTrack(:,2).*obj.pixelLength;
    y = thisTrack(:,3).*obj.pixelLength;
    plot(x,y,'color',color(ii,:)); hold on;
end

set(gca,'fontsize',14,'linewidth',1.5);
xlabel('x [µm]');
ylabel('y [µm]');

axis equal;
pbaspect([1 1 1]);
end
