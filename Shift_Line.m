function shiftedLines = Shift_Line(centerLine, shiftTargets, inds2sample)
%{
-About-
shifts a cell centerlin (or any N-by-2 series of linesegments) orthogonal
to its axis by distances specified by shiftTargets.

-Inputs-
centerLine: either the centerLine of a cell or a cell mesh. if a cell
       mesh is provided, it will be converted into a center line

shiftTargets: points at which the centerLine should be resampled at

inds2sample: indices of the centerLine that should be resampled.  if not
       provided, all indices will be sampled according to input "shiftTargets"

-Outputs-
The shifted result(s)


-Author-
Brad Parry
%}


% shiftTargets = [-3,-2,-1,0,1,2,3]; %distances to re-sample the centerline at

if nargin == 2
    inds2sample = 1:size(centerLine,1);
end
inds2sample = inds2sample(:)';

%Rotation matrix to get an orthogonal line segment
R = [0,-1;1,0];

%initialize the output container
shiftedLines = cell(1,length(shiftTargets));
shiftedLines = cellfun(@(x) zeros(size(centerLine))+NaN, shiftedLines, 'uniformOutput',false);

for k = inds2sample  
    
    if k == 1
        %if the index is one, cant span adjacent center line points,
        %instead only use the first and second points
        u = [centerLine(k,:); centerLine(k+1,:)];
    elseif k == size(centerLine,1)
        %if the index is the final index in the center line, cant span
        %adjacent center line points, instead use the last and second to
        %last
        u = [centerLine(k-1,:); centerLine(k,:)];
    else
        u = [centerLine(k-1,:); centerLine(k+1,:)];
    end
    
    %Sorry for the rough nomencenterLineature here...the variables in caps are
    %meant to be ~static during the loop. Their lower-case relatives can
    %change on each iteration. M represents slope.
    M = diff(u*R);
    MNORM = norm(M);
    
    sampleIndex = 1;
    for dx = shiftTargets
        %scale by the indicated shift
        m = M / MNORM*dx;

        shiftedLines{sampleIndex}(k,:) = centerLine(k,:) + m;
        sampleIndex = sampleIndex + 1;
    end
    
end

shiftedLines = cellfun(@(x) x(~isnan( sum(x,2) ),:), shiftedLines, 'uniformoutput',false);
end