%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% function function data = findPeakLocation(cellList,distance,minPeakHeight,disp,signal)
% @author: Ahmad Paintdakhil and Molly Scott
% @date: November 6, 2015
%=========================================================================
% ********************** input **********************
%cellList: CellList output from Oufti, probably from individual frames 
%distance: minimum distance (in pixels) between peaks of signal intensity
%to be called a peak by the function
%minPeakHeight: Threshold of fluorescence intensity to be called a peak
%disp: Set to 1 to see a plot of the intensity profile for all cells,
%normalized by cell length, where the peaks will be labeled by a triangle.
%Can be quality control.
%signal: Channel for the fluorescence, to be input as [1,0]-->signal 1 or [0,1]-->signal 2 
% ********************** output **********************
%data = cell structure for each cell in the dataset that contains the
%fields pks (peak fluorescence values), locs (locations of the peaks, in
%pixels), normCellLength (normalized cell length, for plotting), signal
%(values of the fluorescence intensity across the whole cell length), and
%cellLength (in pixels)
%=========================================================================
%This script serves to output the peaks of fluorescence intensity in a
%dataset, as well as their locations along the cell length. It allows you
%to plot the data for all cells, in order to have some quality control for
%determining your fluorescence threshold.
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

function data = findPeakLocation(cellList,distance,minPeakHeight,disp,signal)

%data = findPeakLocations(cellList,20,1.3,1,'signal2');
allCells = cat(1,cellList.meshData);
allCells = cat(2,allCells{:});
data = cell(1,length(allCells));
for ii = 1:length(allCells)
    if ~isempty(allCells{ii})
        if signal(1) == 1 && sum(signal) == 1;
            if isfield(allCells{ii},'signal1') && ~isempty(allCells{ii}.signal1)
                signalValue = allCells{ii}.signal1;
            end
        elseif signal(2) == 1 && sum(signal) == 1;
            if isfield(allCells{ii},'signal2') && ~isempty(allCells{ii}.signal2)
                signalValue = allCells{ii}.signal2;
            end
        
        end
        %signalValue = (signalValue - min(signalValue))/(max(signalValue) - min(signalValue));
        signalValue = signalValue/mean(signalValue);
        signalValue = smooth(signalValue,5);
        try
            [pks,locs] = findpeaks(double(signalValue),'minpeakheight',minPeakHeight,'MINPEAKDISTANCE',min(distance,length(signalValue)));
        catch
            continue;
        end
        if isfield(allCells{ii},'lengthvector')
            cellLength = allCells{ii}.lengthvector;
            normCellLength = (cellLength - min(cellLength))/(max(cellLength) - min(cellLength));
            data{ii}.pks = pks;
            data{ii}.locs = locs;
            data{ii}.normCellLength = normCellLength;
            data{ii}.signal = signalValue;
            data{ii}.cellLength = allCells{ii}.length;
        end
    end
end
if disp == 1
    figure;
	hold on
    xx = jet(length(data));
        for ii = 1:length(data)
            try
                plot(data{ii}.normCellLength,abs(data{ii}.signal),'color',[xx(ii,1) xx(ii,2) xx(ii,3)]);%plot curves
				
                plot(data{ii}.normCellLength(data{ii}.locs),data{ii}.pks,'k^','markerfacecolor',[xx(ii,1) xx(ii,2) xx(ii,3)]);%plot peaks
            catch 
            continue;
            end
        end

end
set(gca,'FontSize',14)
xlabel('Relative cell length')
ylabel('Fluorescence (a.u.)')
%set(gca,'Ylim',[0 10]);