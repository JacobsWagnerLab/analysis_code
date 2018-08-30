%% Fit MSD curves based on the specific methods (linear, loglog, ensemble).
%
% -About-
%    Fit the MSD curves (individual or ensemble-averaged) using various
%    methods. The results from the fitting are then used to calculate the
%    diffusion coefficients and alpha values depending on the fitting
%    method used.
%
% -Inputs-
%   - obj:      spt object
%   - fitRange: range of fit (e.g. 1:10, first ten points)
%   - method:   linear, loglog, ensemble, specified as a string
%
% -Output-
%   - D:      diffusion coefficients [µm^2/s]
%   - alpha:  alpha values
%   - alpha:  (ensemble-averaged MSD curve) Only prints the result.
%
% -Example-
%   % Calculate the diffusion coefficient by fitting the linear MSDs
%   - myParticle.fitMSD(1:10,'linear');
%   % Calculate the alpha values by fitting the loglog MSDs
%   - myParticle.fitMSD(2:11,'loglog');
%   % Print out the alpha values by fitting the ensemble loglog MSD
%   - myParticle.fitMSD(2:11,'ensemble');
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function fitMSD(obj,fitRange,method)

% If no fitRange is specified, the range is determined by the shortest
% track. If no mode is specified, the default is linear.
if nargin < 2
    % Calculate the track length
    obj.getTrackLen();
    % Set the default fit range to be the shortest track length
    fitRange = 1:min(obj.trackLen);
    % Default fitting method
    method = 'linear';
elseif nargin < 3
    % Default fitting method
    method = 'linear';
end

% If no MSD or ensemble-averaged MSD were calculated and stored before,
% calculate them now.
if isempty(obj.MSD) || isempty(obj.meanMSD)
    % Calculate the MSD
    obj.getMSD();
    % Calculate the ensemble-averaged MSD
    obj.getMeanMSD();
end

% Switching between different fitting methods
switch lower(method)
    
    % Fit linear MSD to determine diffusion coefficients
    case 'linear'
        % Memory for diffusion coefficients
        D = zeros(length(obj.MSD),1);
        % Loop through MSD curves
        for ii = 1 : length(obj.MSD)
            % Extract a single MSD curve
            thisMSD = obj.MSD{ii};
            % Convert the time and MSD into physical units
            t = thisMSD(:,1).*obj.frameTime;
            msd = thisMSD(:,2).*obj.pixelLength^2;
            % Linear regression based on the fit range
            fitRes = polyfit(t(fitRange),msd(fitRange),1);
            % Store the quarter of the slope as the diffusion coefficient
            % Assumed 2D particle tracking here
            D(ii) = fitRes(1)/4; % [µm^2/s]
        end
        % Pass the result into the instance
        obj.D = D;
    
    % Fit loglog MSD to determine alpha values
    case 'loglog'
        % Memory for alpha values
        alpha = zeros(length(obj.MSD),1);
        % Loop through MSD curves
        for ii = 1 : length(obj.MSD)
            % Extract a single MSD curve
            thisMSD = obj.MSD{ii};
            % Convert the time and MSD into physical units
            t = thisMSD(:,1).*obj.frameTime;
            msd = thisMSD(:,2).*obj.pixelLength^2;
            % If the left bound of the fitRange is 1, warn the user that
            % it is mathematically incorrect and then update the left bound
            % to 2.
            if (fitRange(1) == 1)
                fitRange(1) = 2;
                warning('loglog fit can''t operate on delay 0, skip the 1st point now...')
            end
            % Linear regression based on the fit range
            fitRes = polyfit(log10(t(fitRange)),log10(msd(fitRange)),1);
            alpha(ii) = fitRes(1);
        end
        % Pass the result into the instance
        obj.alpha = alpha;
    
    % Fit the loglog ensemble-averaged MSD curve
    case 'ensemble'
        % If the left bound of the fitRange is 1, warn the user that
        % it is mathematically incorrect and then update the left bound
        % to 2.
        if (fitRange(1) == 1)
            fitRange(1) = 2;
            warning('loglog fit can''t operate on delay 0, skip the 1st point now...')
        end
        
        % Convert the time and MSD into physical units
        x = obj.meanMSD(:,1).*obj.frameTime;
        y = obj.meanMSD(:,2).*obj.pixelLength^2;
        
        % Linear regression
        fitRes = polyfit(log10(x(fitRange)),log10(y(fitRange)),1);
        res = fitRes(1);
        
        % Print the result to the command window
        fprintf('The alpha of the ensemble-averaged MSD curve is %.2f.\n',res);
    
    % If any other key was used for specifying the fitting methods, raise the error.     
    otherwise
        error('Invalid fitting method. Valid options: linear|loglog|ensemble.');
end

end