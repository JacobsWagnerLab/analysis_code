function [] = PlotDotDensity(raw_x,raw_y,xAxis_label,yAxis_label,radius)
%{
-About-
This script will draw plot of x,y 2D density dot plot overlaid with binned mean/ste in red.
Usage example:
PlotDotDensity(raw_x,raw_y,'RNAP travel time (sec)','Ribosome number',0.5)

-Inputs-
(raw_x,raw_y) raw data for (x,y)
xAxis_label : label for the x axis in the plot
yAxis_label : label for the y axis in the plot
radius = default 0.5

-varargin-

-Outputs-
plot of x,y 2D density dot plot overlaid with binned mean/ste in red

-Example-
Format_Example(randn(10,10),randn(10),'error_lim',1)
   
-Supplementary-
Supplemental file location

-Keywords-
2d density plot, bin data

-Dependencies-

-References-
Example plot = Fig. 1F in Campos, et al. Cell (2014)

-Author-
Manuel Campos (original) Sangjin Kim (modified and submitted), 2018 August 29
%}

% If there are too many (raw_x, raw_y), the plot becomes too large in size.
% choose only 4000 pairs randomly but one can change 4000 value to
% something else.
Xraw = [raw_x, raw_y];
tmp = unidrnd(size(Xraw,1),1,4000); 
Xtmp = Xraw(tmp,:);

% this section is from 'denScatter3' from Manuel Campos
h = denScatter3(Xtmp(:,1),Xtmp(:,2),0.5,xAxis_label,yAxis_label);

% See bin_fixBin in the shareRepo
% Alternatively, 
% [binxMean, binyMean,binN, stdy]=bin_fixBin(0,Xraw(:,1),Xraw(:,2),[],[],max(Xraw(:,1)),min(Xraw(:,1)),sqrt(length(Xraw(:,1))));
% One can change the last three input values for bin_fixBin depending on
% the purpose. See bin_fixBin for other possibilities
[binxMean, binyMean,binN, stdy]=bin_fixBin(0,Xtmp(:,1),Xtmp(:,2),[],[],max(Xtmp(:,1)),min(Xtmp(:,1)),sqrt(length(Xtmp(:,1))));

xname={xAxis_label};
yname={yAxis_label};
X = Xtmp(:,1);
Y = Xtmp(:,2);
normX=(X-min(X))./(max(X)-min(X));
normY=(Y-min(Y))./(max(Y)-min(Y));

% this is to determine radius of neighborhood to calculate density
Dlin=pdist([normX,normY],'euclidean');Dsq=squareform(Dlin);
sD=sort(Dlin);
radiusD=sD(round(radius*length(sD)));
neighD=sum(Dsq>=radiusD,2);
neighD=neighD-min(neighD)+1;
colorm = repmat((0.5:.4/max(neighD):1)', 1, 3);
colors=colorm(neighD,:);
% For a better visual effect, plot the denser color at the end:
[~,IXn]=sort(neighD,'descend');
X=X(IXn);Y=Y(IXn);colors=colors(IXn,:);

% first plot shows scatter version of the raw data only
h0 = figure;
scatter(X,Y,15,colors,'filled');
set(gca,'fontname','arial','fontsize',16,'xcolor','k','ycolor','k');
xlabel(xname,'fontname','arial','fontsize',16,'color','k');
ylabel(yname,'fontname','arial','fontsize',16,'color','k');
colormap(flipud(colorm));
hlg=colorbar;
clbLim = get(hlg,'Limits');
set(hlg,'YTick',clbLim,'YtIckLabel',{'1', num2str(max(neighD))},'fontname',...
'arial','fontsize',14,'xcolor','k','ycolor','k');

% second plot shows scatter version of the raw data + binned mean/ste
% together
h = figure;
scatter(X,Y,15,colors,'filled');
hold on;
errorbar(binxMean,binyMean,stdy,'o','linewidth',1,'color',[1 .2 0],'markerfacecolor',[1 .2 0],'markersize',6);% 
set(gca,'fontname','arial','fontsize',16,'xcolor','k','ycolor','k');
xlabel(xname,'fontname','arial','fontsize',16,'color','k');
ylabel(yname,'fontname','arial','fontsize',16,'color','k');
colormap(flipud(colorm));
hlg=colorbar;
clbLim = get(hlg,'Limits');
set(hlg,'YTick',clbLim,'YtIckLabel',{'1', num2str(max(neighD))},'fontname',...
'arial','fontsize',14,'xcolor','k','ycolor','k');
end

function [binxMean, binyMean,binMemnum,SEM]=bin_fixBin(flag,x,y,x2,y2,tE,bE,nB)
% x and y for input data with colors
% tE = topEdge
% bE = botDdge
% nB = numbins
% flag:     if flag is 1, process only first datasets. Otherwise process two
% datasets.
if flag ~= 1
    binEdges = linspace(bE, tE, nB+1);
    [h,whichBin] = histc(x, binEdges);% h give me the number of samples which are in bin nuber i.
    % h(end) is the sample number which x value == binEdges(end);

    flagBinMembers=[];%If their bin is i, show 1, otherwise 0.Logical vecter.
    binMembers=[];% Extracted Y values which X is in this bin.
    binxMean=[];% mean of the binned x value
    binyMean=[];% mean of y value which is corresponding to binnd X group
    binMemnum = [];% number of members in that bin
    SEM = [];
    for i = 1:nB
        flagBinMembers = (whichBin == i);%If their bin is i, show 1, otherwise 0.Logical vecter.
        binMembers     = y(flagBinMembers);% Extracted Y values which X is in this bin.
        binMemnum(i) = length(binMembers);
        binyMean(i)     = nanmean(binMembers);% mean of y value which is corresponding to binnd X group
        binxMean(i)  = mean([binEdges(i) binEdges(i+1)]);% mean of the binned x value
        SEM(i) = nanstd(binMembers);%/(length(binMembers))^0.5;
    end
    for l = length(binMemnum):-1:1
        if binMemnum(l)<=4
             binxMean(l)=[];
             binyMean(l)=[];
             SEM(l)=[];
        end
    end

    binEdges = linspace(bE, tE, nB+1);
    [h,whichBin] = histc(x2, binEdges);% h give me the number of samples which are in bin nuber i.
    % h(end) is the sample number which x value == binEdges(end);

    flagBinMembers=[];%If their bin is i, show 1, otherwise 0.Logical vecter.
    binMembers=[];% Extracted Y values which X is in this bin.
    binxMean2=[];% mean of the binned x value
    binyMean2=[];% mean of y value which is corresponding to binnd X group
    binMemnum = [];% number of members in that bin
    stdy2 = [];
    for i = 1:nB
        flagBinMembers = (whichBin == i);%If their bin is i, show 1, otherwise 0.Logical vecter.
        binMembers     = y2(flagBinMembers);% Extracted Y values which X is in this bin.
        binMemnum(i) = length(binMembers);
        binyMean2(i)     = nanmean(binMembers);% mean of y value which is corresponding to binnd X group
        binxMean2(i)  = mean([binEdges(i) binEdges(i+1)]);% mean of the binned x value
        stdy2(i) = nanstd(binMembers);
    end
    for l = length(binMemnum):-1:1
        if binMemnum(l)<=4
             binxMean2(l)=[];
             binyMean2(l)=[];
             stdy2(l)=[];
        end
    end
    
    axes1 = axes('Parent',figure);
    hold(axes1,'all');
    hold on
    hold on
    plot(x2,y2,'MarkerSize',15,'Marker','.','LineStyle','none','Color',[0.70 0.78 1])
    plot(x,y,'MarkerSize',15,'Marker','.','LineStyle','none','Color',[0.5 0.5 0.5])
    plot(binxMean, binyMean,'MarkerSize',30,'Marker','.','LineStyle','none','Color',[0.2 0.2 0.2])
    errorbar(binxMean, binyMean,SEM,'LineStyle','none','Color',[0.2 0.2 0.2])
    plot(binxMean2, binyMean2,'MarkerSize',10,'Marker','o','LineStyle','none','Color',[0 0 1])
    errorbar(binxMean2, binyMean2,stdy2,'LineStyle','none','Color',[0 0 1])
else
    binEdges = linspace(bE, tE, nB+1);
    [h,whichBin] = histc(x, binEdges);% h give me the number of samples which are in bin nuber i.
    % h(end) is the sample number which x value == binEdges(end);

    flagBinMembers=[];%If their bin is i, show 1, otherwise 0.Logical vecter.
    binMembers=[];% Extracted Y values which X is in this bin.
    binxMean=[];% mean of the binned x value
    binyMean=[];% mean of y value which is corresponding to binnd X group
    binMemnum = [];% number of members in that bin
    SEM = [];
    for i = 1:nB
        flagBinMembers = (whichBin == i);%If their bin is i, show 1, otherwise 0.Logical vecter.
        binMembers     = y(flagBinMembers);% Extracted Y values which X is in this bin.
        binMemnum(i) = length(binMembers);
        binyMean(i)     = nanmean(binMembers);% mean of y value which is corresponding to binnd X group
        binxMean(i)  = mean([binEdges(i) binEdges(i+1)]);% mean of the binned x value
        SEM(i) = nanstd(binMembers);%/(length(binMembers))^0.5;
    end
    for l = length(binMemnum):-1:1
        if binMemnum(l)<=4
             binxMean(l)=[];
             binyMean(l)=[];
             SEM(l)=[];
        end
    end
    axes1 = axes('Parent',figure);
    hold(axes1,'all');
    hold on
     plot(x,y,'MarkerSize',15,'Marker','.','LineStyle','none','Color',[0.5 0.5 0.5])
    plot(binxMean, binyMean,'MarkerSize',30,'Marker','.','LineStyle','none','LineWidth',2,'Color',[0.2 0.2 0.2])
    errorbar(binxMean, binyMean,SEM,'LineStyle','none','Color',[0.2 0.2 0.2])
end
end

function h = denScatter3(X,Y,radius,xAxis_label,yAxis_label)
[binxMean, binyMean,stdy]=bin_fixBin(0,X,Y,[],[],max(X),min(X),floor(sqrt(length(X)))); %close(gcf);%max(X) min(X)
% [binxMean,binyMean,stdy]=bin_fixN(X,Y,75);stdy=stdy./binyMean;

xname={xAxis_label};
yname={yAxis_label};
normX=(X-min(X))./(max(X)-min(X));
normY=(Y-min(Y))./(max(Y)-min(Y));

Dlin=pdist([normX,normY],'euclidean');Dsq=squareform(Dlin);
sD=sort(Dlin);
radiusD=sD(round(radius*length(sD)));
neighD=sum(Dsq>=radiusD,2);
neighD=neighD-min(neighD)+1;

% topper10=round(0.13*max(neighD));
% colorm=bone(max(neighD)+7*topper10);%set
% colorm=colorm(5*topper10:end-2*topper10+1,:);

colorm = repmat((0.5:.4/max(neighD):1)', 1, 3);
colors=colorm(neighD,:);
% For a better visual effect, plot the denser color at the end:
[~,IXn]=sort(neighD,'descend');
X=X(IXn);Y=Y(IXn);colors=colors(IXn,:);%colorm=colorm(IXn,:);

h0 = figure;
scatter(X,Y,15,colors,'filled');
set(gca,'fontname','arial','fontsize',16,'xcolor','k','ycolor','k');
xlabel(xname,'fontname','arial','fontsize',16,'color','k');
ylabel(yname,'fontname','arial','fontsize',16,'color','k');
colormap(flipud(colorm));
hlg=colorbar;
clbLim = get(hlg,'Limits');
set(hlg,'YTick',clbLim,'YtIckLabel',{'1', num2str(max(neighD))},'fontname',...
'arial','fontsize',14,'xcolor','k','ycolor','k');

h = figure;
scatter(X,Y,15,colors,'filled');
% hold on;

hold on;
errorbar(binxMean,binyMean,stdy,'o','linewidth',1,'color',[1 .2 0],'markerfacecolor',[1 .2 0],'markersize',6);% 'Orangee' [1 .4 .2]
set(gca,'fontname','arial','fontsize',16,'xcolor','k','ycolor','k');
xlabel(xname,'fontname','arial','fontsize',16,'color','k');
ylabel(yname,'fontname','arial','fontsize',16,'color','k');
colormap(flipud(colorm));
% h=colorbar('YTick',[1 9.8],'YtIckLabel',{'1', num2str(max(neighD))},'fontname','arial','fontsize',14);%[1 length(colorm)]
% ylabel(h,'density','fontname','arial','fontsize',14,'fontweight','b','color','k');
hlg=colorbar;
clbLim = get(hlg,'Limits');
set(hlg,'YTick',clbLim,'YtIckLabel',{'1', num2str(max(neighD))},'fontname',...
'arial','fontsize',14,'xcolor','k','ycolor','k');
end
