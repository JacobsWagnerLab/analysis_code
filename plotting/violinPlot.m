%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% function varargout = violinPlot(G,data,pix2mu,error_bar,direction)
% @author: Molly Scott (adapted from Manuel Campos)
% @date: November 7, 2015
%=========================================================================
% ********************** input **********************
%G: an array of strings that defines the number of discrete groups in your
%dataset (i.e. 'One zone','Two Zone'....)
%data: a cell array organized by your groups, such that data{1} contains
%your data from your first group, data{2} contains your data for your
%second group, etc.
%pix2mu: pixel to micron conversion factor
%error_bar: if 1, show +/- SEM, if 0, show +/- std deviation
%direction: if 1, will plot by group on y axis, if 0 will plot by group on
%x axis
% ********************** output **********************
%violin plot that shows a scattered density for each group of data, equally
%spaced apart on the y axis.
%=========================================================================
%This script allows you to scatter plot points such that you can visualize
%density of data. It separates discrete scatter plots on the same plot, as
%defined by the user (groups).
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

function varargout = violinPlot(G,data,pix2mu,error_bar,direction)


% Define how many groups are present
unNb=unique(G);
maxcol=length(unNb);
% Define color scale
clmap=winter(maxcol);
if maxcol<3 
    clmap=prism(6);
    clmap(1,:)=clmap(4,:);
    clmap(2,:)=clmap(5,:);
end


M=zeros(length(unNb),1);S=zeros(length(unNb),1);
h=figure;hold on;
% Loop through features (Y columns) and loop to plot groups (G) for each
% feature

  
    for ii = 1 : length(unNb)
        gp=linspace(1,length(unNb),length(unNb)); %creating linear space so that the error bars are properly spaced on x axis
        tmpData=(data{ii}); %each #zones dataset
        tmpData=tmpData*pix2mu; %convert to microns
        [f,xi] = ksdensity(tmpData);
        jitter = interp1(xi,f,tmpData);
        jx=0.25*((randn(size(tmpData))).*jitter./(2*max(jitter))); %how much to jitter on x axis so that we can visualize all data
        plot(tmpData,ii+jx,'o','markeredgecolor',clmap(ii,:),...
            'markerfacecolor',clmap(ii,:)+(1-clmap(ii,:))*0.65);
        if error_bar == 1
            M(ii)=nanmean(tmpData);S(ii)=nanstd(tmpData)/sqrt(length(tmpData)); %mean and SEM
        elseif error_bar == 0
            M(ii)=nanmean(tmpData);S(ii)=nanstd(tmpData);
        end
    end
%     errbar(M,gp,S,'r+-','linewidth',3,'horiz');
    if direction == 1
        plot(M,gp,'ok','markersize',8,'markerfacecolor','k');
    elseif direction == 0
        plot(gp,M,'ok','markersize',8,'markerfacecolor','k');
    end


set(gca,'fontsize',16,'fontname','arial');box off;

% ylabel('Number of zones','fontsize',16,'fontname','arial');
% xlabel('Cell length (\mum)','fontsize',16,'fontname','arial');
% axis([0 max(data{length(unNb)})*pix2mu+10 0 length(unNb)+1])

if nargout==1
    varargout{1}=h;
end
end