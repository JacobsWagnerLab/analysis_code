function center_line = Taylor_Smooth(center_line, varargin)
%{
-About-
Taylor_Smooth uses local first and second derivatives to smooth a cell
centerline (although any N-by-2 matrix may be passed) at its end points.

-Inputs-
center_line: an N-by-2 matrix

-varargin-
'adjustment_start_index'  set the limit of acceptable error in the analysis

'adjustment_length'       the number of segments that should be adjusted.
    if ajdustment_length is larger than adjustment_start_index, the new
    centerline will extend beyond the original line

'taylor_region'     the coordinates to be used for derivative estimation in
    the taylor series expansion

-Outputs-
centerline:     the center_line with each end smoothed

-Example-
   
-Supplementary-

-Keywords-
cell mesh pole smooth pole 

-Dependencies-

-References-

-Author-
Brad Parry
%}

adjustment_start_index = 7;
adjustment_length = 7;
taylor_region = 10:-1:2;

for k = 1:length(varargin)
    if strcmpi(varargin{k},'adjustment_start_index')
        adjustment_start_index = varargin{k+1};
    end
    
    if strcmpi(varargin{k},'adjustment_length')
        adjustment_length = varargin{k+1};
    end
    
    if strcmpi(varargin{k},'taylor_region')
        taylor_region = varargin{k+1};
    end
end

region_to_smooth = (adjustment_length-1):-1:1;

tmp0 = adjust_one_side(center_line, taylor_region, adjustment_start_index, region_to_smooth);
center_line = flipud( adjust_one_side( flipud(tmp0), taylor_region, adjustment_start_index, region_to_smooth) );
end

function out = adjust_one_side(input_xy, taylor_region, adjustment_start_index, region_to_smooth)
polar_derivatives = mean(diff(input_xy(taylor_region,:)));
polar_derivatives2 = mean(diff(diff(input_xy(taylor_region,:))));

input_xy = input_xy(adjustment_start_index:end,:);

modification_values=[];
modification_values(:,1) = input_xy(1,1) + region_to_smooth*polar_derivatives(1)' + (region_to_smooth.^2)/2*polar_derivatives2(1);
modification_values(:,2) = input_xy(1,2) + region_to_smooth*polar_derivatives(2)' + (region_to_smooth.^2)/2*polar_derivatives2(2);

out = cat(1, modification_values, input_xy);
end
