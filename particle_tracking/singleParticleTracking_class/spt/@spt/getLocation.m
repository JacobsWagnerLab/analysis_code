%% Determine particle locations.
%
% -About-
%   The first step in particle tracking relies on the accurate
%   localizations of the objects within an image frame. The following
%   procedures are implemented by the method getLocation():
%       1. Filter a frame with a band-pass filter
%       2. Set any pixels with intensities below a pre-defined threshold to
%          zero
%       3. Create a binary copy of the frame and segregate the objects
%       4. Eliminate any objects with sizes outside the pre-defined size
%          range
%       5. Loop through individual objects and create a copy of the object
%          image with extra padding
%       6. Localize the brightest spot within the object
%       7. Crop out a 9x9 image with the brightest spot at the center
%       8. Fit the 9x9 image with a bivariate Gaussian distribution to
%          determine the particle location as Gaussian peak position
%       9. Calculate the actual particle position in reference to the
%          original full frame
%
%   The method has defaults built in for all the required parameters. The
%   user can use keywords to overwrite these default in order to set the
%   localization properly. See examples below.
%
%   The method also assumes the user has acquired images in a 16-bit
%   camera, from with the brightest pixel presented can only be 2^16-1. If
%   any pixel with intensity higher than this value is detected, the entire
%   frame will be skipped and the user would be warned about the detection
%   of a saturated image. A saturated image can result inaccurate fitting
%   with the Gaussian distribution.
% 
% -Input- 
%   - obj: spt object
% 
% -Varargin- 
%   - minSpotSize: the smallest object to be localized (default: 0) 
%   - maxSpotSize: the largest object to be localized (default: Inf)
%   - intensityThreshold: the intensity threshold below which a pixel will
%     be zeroed after the band-pass filter (default: 0)
%   - startFrame: the first frame to be analyzed (default: 1) 
%   - endFrame: the last frame to be analyzed (default: all frames)
% 
% -Output- 
%   - locs: a Nx3 matrix with all particle positions localized. N
%           is the number of particle positions localized, and the three
%           columns represent particle positions in x, y dimensions and
%           frame numbers at which the particles are localized,
%           respectively.
% 
% -Example-
%   % Localize particles with default parameters
%   myParticle.getLocation()
%   % Localize particles with user-defined parameters
%   myParticle.getLocation('minSpotSize',        10, ...
%                          'maxSpotSize',        30, ...
%                          'intensityThreshold', 10, ...
%                          'startFrame',         1, ...
%                          'endFrame',           1000);
%
% -Author- 
%   Yingjie Xiang, CJW Lab, Yale University


function getLocation(obj, varargin)

% Default parameter values
% Default:
%   - localize for the entire stack
%   - localize objects with a size between 1 to inifite pixels
%   - reset only negative pixels after convolution

% Default parameters
% Minimum object size to be localized
minSpotSize = 1;
% Maximum object size to be localized
maxSpotSize = Inf;
% Intensity threshold, below which the pixel will be set to 0
intensityThreshold = 0;
% First frame to be analyzed
startFrame = 1;
% Last frame to be analyzed
endFrame = obj.imNum;

% Parse user inputs and overwrite the defaults
for ii = 1: 2: length(varargin)
    switch lower(varargin{ii})
        case 'minspotsize'
            minSpotSize = varargin{ii+1};
        case 'maxspotsize'
            maxSpotSize = varargin{ii+1};
        case 'intensitythreshold'
            intensityThreshold = varargin{ii+1};
        case 'startframe'
            startFrame = varargin{ii+1};
        case 'endframe'
            endFrame = varargin{ii+1};
    end
end

% Store the localization parameters
obj.iThreshold  = intensityThreshold;
obj.minSpotSize = minSpotSize;
obj.maxSpotSize = maxSpotSize;

% Memory for localized particle positions
locs = [];
for ii = startFrame:endFrame
    
    % Extract a single frame
    thisFrame = double(imread(obj.stackPath,'Index',ii));
    
    % Filter the image using a band-pass filter
    % The intensity threshold is specified, so that after the band-pass,
    % any pixels with intensities lower than the threshold will be zeroed.
    bandPassedIm = bpass(thisFrame,1,obj.estimatedSize,obj.iThreshold);
    % Create binary images for object segregations
    bandPassedIm(bandPassedIm>0) = 1;
    % Get connected components
    bw = bwconncomp(bandPassedIm);
    
    % Remove objects outside the pre-set range (minSpotSize, maxSpotSize)
    % Number of pixels in each object
    sz = cellfun(@length, bw.PixelIdxList);
    % Index the objects with sizes outside the favorable range
    kill = (sz<obj.minSpotSize) | (sz>obj.maxSpotSize);
    % Eliminate the objects outside range
    bw.PixelIdxList(kill) = [];
    % Update the number of objects
    bw.NumObjects = length(bw.PixelIdxList);
    
    % Memory for storing the positions of objects in one frame
    spots = [];
    % The pad is needed for the later 2D Gaussian fitting
    pad = 3;
    % Loop through each object
    for jj = 1:length(bw.PixelIdxList)
        % For each object, find a sub-image 
        [r,c] = ind2sub(bw.ImageSize, bw.PixelIdxList{jj});
        % Create a sub-image of this specific object with an extra of 3
        % pixels of pad
        sub_im = thisFrame(min(r)-pad:pad+max(r),min(c)-pad:pad+max(c));
        % Warn the user if a saturated object is detected and skip this
        % object
        if max(sub_im(:)) == 2^16-1
            warning('Saturated object dectected (16 bit), now skipping...');
            continue
        else    
            % Find the brightest pixel in the sub-image
            [maxR, maxC] = find(sub_im == max(sub_im(:)));
            % Position of the brightest pixel
            maxR = maxR(1);
            maxC = maxC(1);
            
            % If the brightest pixel is very close to the boundary, we
            % cannot create the truncated image described below, since we
            % may need pixels beyond the boundary of the sub-images
            % (non-existing).
            if (maxR-4 < 1) || (maxR+4 > size(sub_im,1)) || ...
                    (maxC-4 < 1) || (maxC+4 > size(sub_im,2))
                continue;
            end
            
            % Create 9x9 trucated image with the brightest pixel at the
            % center
            t_im = sub_im(maxR-4:maxR+4, maxC-4:maxC+4);
            
            % Fit the truacted image with a bivariate Gaussian distribution
            % and output the position of the peak (i.e. particle position).
            % Note: this position is relative to the trucated-image.
            sfit = gaussfit2d(t_im);
            xy = [sfit.x0,sfit.y0];
            
            % Back-calculate the particle position in the refernece frame
            % of entire image
            x = xy(1,1)+ maxC-4 -1 + min(c) - pad - 1;
            y = xy(1,2)+ maxR-4 -1 + min(r) - pad - 1;
            
        end
        
        % Spots holds two columns (positions in x and y dimensions)
        spots(end+1,1:2) = [x,y];
    end
    
    % To account some replications in localized positions
    spots = unique(spots,'rows');
    % Append a column of frame numbers
    loc = [spots, ones(size(spots,1),1).*ii];
    % Append to store all particle locations
    locs = [locs; loc];
end

% Pass the results to the object
obj.locs = locs;


end