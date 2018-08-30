function [interior, exterior] = resizeCellMesh(cellMesh, shiftOffset)
%{
cellMesh: one Oufti cellMesh (e.g., cellList.meshData{1}{1}.mesh)
shiftOffset: distance from mesh vertices in pixels that should be used for the
    gradient calculation. Reasonable values are in the range of 1 - 2. 

Each vertex in a cell mesh is contracted toward the center of the cell
(output: interior) and toward the exterior of the cell (output: exterior).
This is done by defining a line orthogonal to each vertex scaling it by
shiftOffset. This is done in the code below (getNewVertices),
resizeCellMeshes does geometric housekeeping to ensure that a valid cell
mesh is returned.

Brad Parry. Jacobs-Wagner Lab, 2016
%}
xy = cat(1, cellMesh(:,1:2), flipud(cellMesh(:,3:4)) );

Left = cat(1, cellMesh(2,3:4), cellMesh(:,1:2), cellMesh(end-1,3:4));
Right = cat(1, cellMesh(2,1:2), cellMesh(:,3:4), cellMesh(end-1,1:2));

[L_points1, L_points2] = getNewVertices(Left, shiftOffset);
[R_points1, R_points2] = getNewVertices(Right, shiftOffset);

if ispolycw(xy(:,1), xy(:,2))
    interior = [L_points1, R_points2];
    exterior = [L_points2, R_points1];
else
    interior = [L_points2, R_points1];
    exterior = [L_points1, R_points2];
end

end

function [points1, points2] = getNewVertices(xy, shiftOffset)
%{
linear algebra to contract and expand a polygon from its boundary

Brad Parry. Jacobs-Wagner Lab, 2016
%}
R = [0,-1;1,0];
points1 = zeros(size(xy,1)-2,2);
points2 = points1;
for k = 2:(size(xy,1)-1)
    M = diff( xy([k-1,k+1],:)*R );
    shift1 = M / norm(M)*shiftOffset;
    shift2 = M / norm(M)*(-shiftOffset);
    points1(k-1,:) = xy(k,:) + shift1;
    points2(k-1,:) = xy(k,:) + shift2;
end
end