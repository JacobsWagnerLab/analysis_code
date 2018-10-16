function grad = cellEdgeGradient(cellMesh, shiftOffset, image)
%{
INPUTS
cellMesh: one Oufti cellMesh (e.g., cellList.meshData{1}{1}.mesh)
shiftOffset: distance from mesh vertices in pixels that should be used for the
    gradient calculation. Reasonable values are in the range of 1 - 2. 
image: the phase contrast image that the cellMesh was calculated from

cellEdgeGradient evaluates the image gradient over each vertex in cellMesh

OUTPUT
The gradient (grad) is calculated as an average of the difference in
normalized image values sampled from a dilation and contraction of the cell
mesh.

Brad Parry. Jacobs-Wagner Lab, 2016
%}

%compute a contracted (interior) and dilated (exterior) version of the cell
%mesh with an identical number of vertices
[interior, exterior] = resizeCellMesh(cellMesh, shiftOffset);
%convert the new meshes into a polygon
interior = cat(1, interior(:,1:2), flipud(interior(:,3:4)) );
exterior = cat(1, exterior(:,1:2), flipud(exterior(:,3:4)) );
%round the polygons so that they can be used for indexing.
%***
%***    it may be more interesting/accurate to sample the image by
%***    interpolation. to do so, the next lines would be replaced with an
%***    interpolation routine.
%***
interior = round(interior);
exterior = round(exterior);
try
    interior = sub2ind(size(image), interior(:,2), interior(:,1));
    exterior = sub2ind(size(image), exterior(:,2), exterior(:,1));
    grad = mean( double(image(exterior)) - double(image(interior)) );
catch
    grad = NaN;
end
end