function iLine = Project_Cell_Signal(cell_mesh, shift_targets, pole_length, image, varargin)
%{
-About-
This function projects image intensity values into cellular coordinates by
projecting pixel values onto the cell centerline. The cell centerline is
calculated internally from an Oufti determined cell mesh
(i.e., cellList.meshData{frame}{cell}.mesh).

Project_Cell_Signal will calculate the cell centerline, and shift as many
copies of it as necessary to satisfy shift_targets. These shifted
centerlines are then used to sample the input image and an output matrix of
image intensities under the shifted centerlines will be output in a matrix
of size [cell length - 2*pole_length, length(shift_targets)]

-Inputs-
cell_mesh: an Oufti cell mesh; cell_mesh = cellList.meshData{F}{C}.mesh;

shift_targets: indicates pixel deviations from the cell's centerline that
should be sampled shift_targets = -2:2; Shifts are performed in each 
direction orthogonal to the centerline.

pole_length: The number of pixels should be excluded from each cell pole in
the analysis

image: is the raw image that the cell_mesh corresponds to

-varargin-
center_line_only: if an N-by-2 polygon is passed as a cell mesh,
center_line_only should be set to True. In this case, the image pixels will
be projected onto the polygon provided as cell_mesh.

Other arguments, if provided, are passed to Taylor_Smooth. See
Taylor_Smooth documentation for additional details on possible varargin

-Outputs-
iLine: will be an N-by-length(shift_targets) matrix of image
intensities projected onto the cell centerline and shifted
representations of the cell centerline as determined by shift_targets

-Example-
   
-Supplementary-

-Keywords-
Cell projections

-Dependencies-
Extract_Cell_Image, Shift_Line, Taylor_Smooth

-References-

-Author-
Brad Parry
%}

center_line_only = false;

for k = 1:length(varargin)
    if strcmpi(varargin{k},'center_line_only')
        center_line_only = true;
        center_line = cell_mesh;
        cellImage=image;
        center_line = Taylor_Smooth(center_line,varargin{:});
    end
end

if center_line_only
    inds2sample = (pole_length+1) : (length(center_line)-pole_length);
else
    [cellImage, bx] = Extract_Cell_Image(cell_mesh, image, 'pad', 10);
    cell_mesh(:,[1,3]) = cell_mesh(:,[1,3]) - bx(3) + 1;
    cell_mesh(:,[2,4]) = cell_mesh(:,[2,4]) - bx(1) + 1;
    
    center_line = cat(2, mean(cell_mesh(:,[1,3]),2), mean(cell_mesh(:,[2,4]),2));
    center_line = Taylor_Smooth(center_line,varargin{:});
    
    inds2sample = (pole_length+1) : (length(center_line)-pole_length);
end

shifted_lines = Shift_Line(center_line, shift_targets, inds2sample);

%set convention that inner curvature will be on the left, outer on the
%right
totalLength = @(x) sum(sum(diff(x).^2,2).^(1/2));
if totalLength(shifted_lines{1}) > totalLength(shifted_lines{end})
    shifted_lines = shifted_lines(length(shifted_lines):-1:1);
end

[y,x] = ndgrid(1:size(cellImage,1),1:size(cellImage,2));
iLine = zeros(size(shifted_lines{1},1), length(shift_targets));
for C = 1:length(shifted_lines)
    x1 = shifted_lines{C}(:,1);
    y1 = shifted_lines{C}(:,2);
    iLine(:,C) = interp2(x,y,double(cellImage),x1,y1,'cubic');
end

end