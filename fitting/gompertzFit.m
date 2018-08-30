function gmpPar = gompertzFit(t, OD, showFit)
%% -------------------------------------------------------------------------
% Gompertz function fitting
% function gmpPar = gompertzFit(t, OD, showFit)
% @auther: Manuel Campos
% @date: August 25 2016
% @copyright 2015-2016 Yale University
%=========================================================================
% ********************** input **********************
%t:         time vector
%OD:        OD data matrix
%showFit:   logical flag indicating whether to show or not each fit. If set
%           to TRUE, hit any key stroke to move on.
%
% ********************** output **********************
%gmpPar     Array of fit parameters. As many lines as fits by 4 columns.
%           [Sat. density; max growth rate; lag time; No]
%
%=========================================================================
% This function is meant to fit a complete growth curve and to extract 3
% parameters from the growth curve and crosscheck with the forth. The
% saturation density, max growth rate and lag time are the first three
% parameters we are seeking. The forth one is the projected OD (below level
% of detection) of the culture at the start of the growth period. A good
% verification would be to check the real OD at saturation (need dilution)
% and verify that this value corrected for the initial dilution matches the
% extimation of the fit
%
% DEVELOPMENT: Return 95% confidence intervals of the parameters
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% Robustness of code vis a vis inputs
if nargin<3 || ~any(ismember(showFit,[0,1]))
    showFit = false;
end
if size(t,1)==1
    t = t';
end
if size(OD,2)==size(t,1)
    OD = OD';
end
if ~isequal(size(t,1),size(OD,1))
    disp('The number of lines in the time vector and the OD data should match');
    gmpPar = [];
    return
end

%% Define the fitting model
ft = fittype( 'a*exp(-exp(b*exp(1)*(c-x)/a + 1)) + d', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Robust = 'Bisquare';
% opts.MaxFunEvals = 1000;


%% Loop through the columns of OD data and fit
gmpPar = zeros(size(OD, 2), 4);
w = waitbar(0, 'Fitting curve 1');
for ii=1:size(OD, 2)
    waitbar(ii/size(OD, 2), w, sprintf('Fitting curve %d',ii));
    % Log transform the normalized data
    grlog = log(OD(:,ii) - 0.087);
    % prepare logicals to remove useless values from fitted data
    ix = ~isinf(grlog) & ~isnan(grlog) & grlog> log(0.09-0.087); 
    if sum(ix)>100
        % Initialize fitting parameters
        init_tlag = t(find(grlog<=-5.5,1,'last'));
        if isempty(init_tlag)
            init_tlag = 50;
        end
        % All parameters (Except teh log of the inital culture density)
        % have to be positive.
        opts.Lower = [0 0 init_tlag/2 -inf];
        opts.StartPoint = [6 0 init_tlag -8];
        
        % Actual fitting and record parameters
        [gmpF, ~] = fit( t(ix), grlog(ix), ft, opts );
        gmpPar(ii,:) = [gmpF.a, gmpF.b, gmpF.c, exp(gmpF.d)];
        % Display if showFit flag set to TRUE
        yfit = feval(gmpF, t(ix));
        if showFit
            figure(2);hold off;
            set(gcf,'position',[50 100 560 420]);
            plot(t, OD(:,ii) - 0.087, '.');hold on;
            plot(t(ix), exp(yfit), '-r');
            set(gca,'fontsize',16,'xcolor','k','ycolor','k','layer','top','yscale','log');
            xlabel('time','fontsize',16,'fontweight','b','color','k');
            ylabel('log OD_{600nm}','fontsize',16,'fontweight','b','color','k');
            title(['Fit curve ',num2str(ii)]);legend off
            pause;
        end
    else
        gmpPar(ii,:) = nan(1, 4);
    end
end
close(w);

end
