function [cellImage, bx] = Extract_Cell_Image(cellMesh, image, varargin)
%{
-About-
returns and image and bounding box around an Oufti cell mesh

-Inputs-
cellMesh: an Oufti cell mesh; cell_mesh = cellList.meshData{F}{C}.mesh;

image: is the raw image that the cell_mesh corresponds to

-varargin-
pad: the number of pixels to pad around the cell boundary

region:

-Outputs-
cellImage: an image cropped around the cell mesh

bx: the bounding box of cellImage in the original image as LBWH

-Author-
Brad Parry
%}


pad = 10;
bx = [];

for k = 1:length(varargin)
    if strcmpi(varargin{k},'pad')
        pad = varargin{k+1};
    elseif strcmpi(varargin{k}, 'region')
        bx = varargin{k+1};
    end
end

cellMesh = double( cat(1, cellMesh(:,1:2), flipud(cellMesh(:,3:4))) );

if isempty(bx)
    bxlim = [1 size(image,1) 1 size(image,2)];
    bx = [floor(min(cellMesh(:,2))), ceil(max(cellMesh(:,2))), floor(min(cellMesh(:,1))), ceil(max(cellMesh(:,1)))];

    bx([1 3]) = bx([1 3]) - pad; 
    bx([2 4]) = bx([2 4]) + pad;
    %if the box is outside of the image, bring it back
    ck = (bx - bxlim).*[1 -1 1 -1];
    bx(ck < 0) = bxlim(ck < 0);
end

%get the cell delimited image
cellImage = image(bx(1):bx(2),bx(3):bx(4));
end