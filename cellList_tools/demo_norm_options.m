%demograph with different signal normalization routines
%-------------------------------------------------------------------------%
%INPUTS: 
%cellList: output from Oufti
%signal: which signal to use
%   -'signal1'
%   -'signal2'
%maxCellLength: estimate of maximum cell length (to curate your dataset)
%normalization: how to initially get signal values for each step length
%along the cell
%   -'area': divide signal by steparea
%   -'volume': divide by stepvolume
%   -'summed_int': divide by summed intensity in the cell
%relative: how to normalize your signal values for display purposes
%   -'by_brightest_perCell': normalize to the maximum fluorescence value in each
%   cell
%   -'by_brightest_pixel': normalize to the maximum fluorescence value in the
%   entire dataset
%   -'by_mean_perCell': normalize, by cell, to that cell's mean
%   fluorescence value
%   -'by_mean_perCell_rel': normalize, by cell, to that cell's mean
%   fluorescence value, then scale from 0 to 1
%   -'by_median_perCell': normalize, by cell, to that cell's median
%   fluorescence value
%   -'by_median_perCell_rel': normalize, by cell, to that cell's median
%   fluorescence value, then scale from 0 to 1
%cmap: colormap
%-------------------------------------------------------------------------%
%OUTPUTS: two plots, (1) demograph oriented so all cells touch the left
%side of the figure and (2) demograph organized at the midcell

function demo_norm_options(cellList,signal,maxCellLength,normalization,relative,cmap)

count = 0;

%count number of cells in dataset
for frames = 1:length(cellList.meshData)
    for cells = 1:length(cellList.meshData{frames})
         if isempty(cellList.meshData{frames}{cells}) || length(cellList.meshData{frames}{cells}.mesh)<4 ...
                    ||~isfield(cellList.meshData{frames}{cells},signal) ...
                    || isempty(cellList.meshData{frames}{cells}.(signal)) || cellList.meshData{frames}{cells}.length>maxCellLength
            continue
        end
            count = count + 1;
    end
end

%create an array with all of your cell lengths, and a structure with all of
%the corresponding signal information
all_lengths = zeros(count,1);
normalized_signals = cell(count,1);
count = 0;
for frames = 1:length(cellList.meshData)
    for cells = 1:length(cellList.meshData{frames})
        if isempty(cellList.meshData{frames}{cells}) || length(cellList.meshData{frames}{cells}.mesh)<4 ...
                    ||~isfield(cellList.meshData{frames}{cells},signal) ...
                    || isempty(cellList.meshData{frames}{cells}.(signal)) || cellList.meshData{frames}{cells}.length>maxCellLength
            continue
        end
            count = count+1;
            all_lengths(count,1) = cellList.meshData{frames}{cells}.length;
            %normalize by area or volume, depending on input
            if strcmp(normalization,'area')
                sgn = cellList.meshData{frames}{cells}.(signal)./cellList.meshData{frames}{cells}.steparea;
            elseif strcmp(normalization,'volume')
                sgn = cellList.meshData{frames}{cells}.(signal)./cellList.meshData{frames}{cells}.stepvolume;
            elseif strcmp(normalization,'summed_int')
                sum_int = cellList.meshData{frames}{cells}.(signal)./cellList.meshData{frames}{cells}.steparea;
                sgn = sum_int./sum(sum_int);
            end
            %determine how you want to normalize for display purposes
            if strcmp(relative,'by_brightest_perCell')
                normalized_signals{count} = sgn / max(sgn);
            elseif strcmp(relative,'by_brightest_pixel')
                normalized_signals{count} = sgn; 
                max_sig_temp = 0;
                for signals = 1:length(normalized_signals)
                    old = max_sig_temp;
                    max_sig_temp = max(normalized_signals{signals});
                    if max_sig_temp > old
                        max_sig_temp = max_sig_temp;
                    else 
                        max_sig_temp = old;
                    end
                end
                normalized_signals{count} = sgn/max_sig_temp;  
            elseif strcmp(relative,'by_mean_perCell')
                normalized_signals{count} = sgn/mean(sgn);
            elseif strcmp(relative,'by_mean_perCell_rel')
                normalized_signals{count} = sgn/mean(sgn);
            elseif strcmp(relative,'by_median_perCell')
                normalized_signals{count} = sgn/median(sgn);
            elseif strcmp(relative,'by_median_perCell_rel')
                normalized_signals{count} = sgn/median(sgn);
            end
    end
end

%if using relative fluor by cell and scaling 0-->1:
if strcmp(relative,'by_mean_perCell_rel') || strcmp(relative,'by_median_perCell_rel')
    for cells = 1:length(normalized_signals)
        normalized_signals{cells} = normalized_signals{cells}/max(normalized_signals{cells});
    end
end

%order your cell length and signal structures by increasing cell length
maxCellLength = max(all_lengths);

s = cellfun(@size,normalized_signals,'uniform',false);
[trash idx] = sortrows(cat(1,s{:}));
ordered_norm_signals = normalized_signals(idx);
all_lengths_ordered = all_lengths(idx);

%create an array with the signal information, by increasing cell length
max_pixels_to_plot = length(ordered_norm_signals{end,1}); %determine the #pixels from the largest cell
cells_to_plot = ones(length(ordered_norm_signals),max_pixels_to_plot); %we want our image matrix to contain ones for the length of the largest cell, and then fill in values for our real data (ones-->background color)

for cells = 1:length(ordered_norm_signals)
    for sig = 1:length(ordered_norm_signals{cells,1})
        cells_to_plot(cells,sig) = ordered_norm_signals{cells,1}(sig,1);
    end
end

%plot cells at left side of image
figure(1)
image(cells_to_plot,'CDataMapping','scaled')
colormap(cmap)

%re-create the array, such that cells are plotted from the midcell
middle = floor(max_pixels_to_plot/2);
cells_to_plot_mid = ones(length(ordered_norm_signals),max_pixels_to_plot);

for cells = 1:length(ordered_norm_signals)-1
    mid_sig = floor(length(ordered_norm_signals{cells,1})/2);
    count = -1;
    for sig = 1:length(ordered_norm_signals{cells,1})
        count = count + 1;
        cells_to_plot_mid(cells,middle-mid_sig+count) = ordered_norm_signals{cells,1}(sig,1);
    end
end

cells_to_plot_mid(length(ordered_norm_signals),:) = ordered_norm_signals{end,1}'; %fill in data for longest cell

%plot cells, centered at the midcell of the longest cell
figure(2)
image(cells_to_plot_mid,'CDataMapping','scaled')
colormap(cmap)

end

    