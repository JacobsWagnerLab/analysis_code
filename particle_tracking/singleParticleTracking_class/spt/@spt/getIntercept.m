%% Find the y-intercepts of all MSD curves.
%
% -About-
%   The diffusion coefficient may also be deduced from the y-intercept of
%   the loglog MSD curves, since MSD=4Dt (2D tracking), log(MSD) = log(4D)
%   + log(t). D = 10^(intercepts)/4. However, in real experiments, we never
%   measure the intercepts (requires infinitely short time delay). Instead,
%   we can either fit the loglog MSD curves using data at multiple delays
%   or simply assume the diffusion is Brownian and extrapolate just using
%   the point at the shortest delay. Here, we implement the latter. The
%   intercepts are deduced in assumption of the normal diffusion at small
%   delay. This assumption is often true for all practical purposes.
%
% -Inputs-
%   - obj: spt object
%
% -Output-
%	- Intercepts: the intercepts of the MSD curves. (Note: the intercepts
%	are in linear scales. In order to calculate the diffusion
%	coefficient, we take log10 transformation.) D = 10^(intercepts) / 4.
%
% -Example-
%   % Calculate the y-intercepts of the loglog MSD
%   myParticle.getIntercept();
%   % Deduce the diffusion coefficients from the intercepts
%   D = 10.^(myParticle.intercepts)./4;
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function getIntercept(obj)
% Memory for collecting the MSD at the smallest delay
x = zeros(length(obj.MSD),1);
% Loop through all MSD curves
for ii = 1 : length(obj.MSD)
    % Extract the MSD at the smallest delay
    x(ii) = obj.MSD{ii}(2,2).*obj.pixelLength^2;
end
% Assuming normal diffusion (loglog MSD slope 1), extrapolate the
% y-intercepts and pass to the object
obj.intercepts = log10(x) - log10(obj.frameTime) * 1;

end