function [x2, y2] = find_inside_boundary(x, y)
% Find points inside the given boundary.
% 
% [x2, y2] = find_inside_boundary(x, y)
% 
% DESCRIPTION:
%   Find the coordinates of points which are normal to the open boundary
%   described by the coordinates (x, y). The distance from the boundary is
%   determined from the mean of the length of the two adjacent boundary
%   element lengths.
% 
% INPUT:
%   x, y - coordinate pairs for the open boundary.
% 
% OUTPUT:
%   x2, y2 - coordinate pairs for the points as normal to the open boundary
%               as possible (i.e. bisecting the angle between the two
%               adjacent open boundary element faces).
% 
% EXAMPLE USAGE:
%   [x2, y2] = find_inside_boundary(x, y)
% 
% NOTES:
%   This works best with cartesian coordinates but will work with spherical
%   too, although the angles for large elements will be incorrect (in an
%   absolute sense).
% 
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
% 
% Revision history:
%   2013-03-11 First version.
% 
%==========================================================================

subname = 'find_inside_boundary';

global ftbverbose
if ftbverbose
    fprintf('\n'); fprintf(['begin : ' subname '\n']);
end

% Check the inputs
if length(x) ~= length(y)
    error('Size of inputs (x, y) do not match.')
else
    % Set the number of points in the boundary
    np = length(x);
end

x2 = nan(np, 1);
y2 = nan(np, 1);

% Order the boundary points in clockwise order (for the polygon centroid
% calculation).
[x1, y1] = poly2cw(x, y);
[cx, cy] = centroid(x1(:), y1(:));

% For each node, find the two closest nodes to it and use those as the
% adjacent vectors. Doing this is slower than just assuming the nodes are
% provided in a sorted order, but means we don't have to worry too much
% about any irregularities from the sorting above (poly2cw).
for pp = 1:np
    
    [~, idx] = sort(sqrt((x1(pp) - x1).^2 + (y1(pp) - y1).^2));
    
    % Get the coordinates of the two nearest points (skip the first closest
    % because that's the current point).
    [px1, py1] = deal(x1(idx(2)), y1(idx(2)));
    [px2, py2] = deal(x1(idx(3)), y1(idx(3)));
    
    % Find the length of the edges of the triangle formed from the three
    % points we're currently considering. 1 and 2 are those adjacent and 3
    % is the remaining side.
    ln1 = sqrt((x1(pp) - px1)^2 + (y1(pp) - py1)^2);
    ln2 = sqrt((x1(pp) - px2)^2 + (y1(pp) - py2)^2);
    ln3 = sqrt((px1 - px2)^2 + (py1 - py2)^2);
    
    % Find the angle between the two element edges and the current node
    % (cosine rule). Use the real component only for cases where the three
    % points lie on a straight line (in which case the angle should be
    % 180 degrees).
    ang1 = real(acosd((ln1^2 + ln2^2 - ln3^2) / (2 * ln1 * ln2)));
    ang1b = ang1 / 2; % bisect the angle
    
    % Find the angle to the horizontal for the current node and one of the
    % other points.
    ang2 = atan2((py1 - y1(pp)), (px1 - x1(pp))) * (180 / pi);
    
    % Find the difference between the two.
    ang3 = ang2 - ang1b;
    
    % Now get the mean length of the two closest element edges and use that
    % to create the new point inside the boundary.
    ml = mean([ln1, ln2]);
    dy = ml * sind(ang3);
    dx = ml * cosd(ang3);

    % Add the offsets to the current node to get the new node's position.
    % This is where things get a bit hairy: we'll assume that the boundary
    % is approximately circular in nature. This means we can use its
    % centroid as a tool to find the points inside the current boundary and
    % those outside.
    [xx(1), yy(1)] = deal(x1(pp) + dx, y1(pp) + dy);
    [xx(2), yy(2)] = deal(x1(pp) - dx, y1(pp) - dy);

    % Find the distances from the centroid.
    dist = sqrt((cx - xx).^2 + (cy - yy).^2);
    
    if dist(1) < dist(2)
        [x2(pp, 1), y2(pp, 1)] = deal(xx(1), yy(1));
%         [x2(pp, 2), y2(pp, 2)] = deal(xx(2), yy(2));
    else
        [x2(pp, 1), y2(pp, 1)] = deal(xx(2), yy(2));
%         [x2(pp, 2), y2(pp, 2)] = deal(xx(1), yy(1));
    end
end

if ftbverbose
    fprintf('end   : %s \n', subname)
end

