function [adjx, adjy] = fix_inside_boundary(x, y, node_ids, ang_thresh)
% Fix unstructured grid points inside the given boundary.
%
% [adjx, adjy] = fix_inside_boundary(x, y, node_ids, ang_thresh)
% 
% DESCRIPTION:
%   Find the coordinates of points which are normal to the open boundary
%   described by the coordinates (x(node_ids), y(node_ids)). The distance
%   from the boundary is determined from the mean of the length of the two
%   adjacent boundary element lengths. Once the 'ideal' position has been
%   identified, find the closest existing nodal position and change it to
%   the 'ideal' position.
%
%   The resulting x2 and y2 coordinates can be exported to 2dm to be
%   checked in SMS for mesh quality with the fvcom-toolbox function
%   write_SMS_2dm.
%
% INPUT:
%   x, y - Unstructured grid coordinates.
%   node_ids - List of IDs of the nodes within the grid which are on the
%              open boundary of interest.
%   ang_thresh - [optional] Specify a minimum angle in degrees between the
%                two adjacent nodal positions to deal with corners better.
%
% OUTPUT:
%   adjx, adjy - New unstructured grid coordinate pairs in which the points
%                just inside the open boundary have been adjusted so as to
%                bisect the angle formed by the two adjacent boundary
%                faces.
%
% EXAMPLE USAGE:
%   [adjx, adjy] = fix_inside_boundary(Mobj.x, Mobj.y, Mobj.read_obc_nodes{1}, 90)
%
% NOTES:
%   This works best with cartesian coordinates but will work with spherical
%   too, although the angles for large elements will be incorrect.
%   Secondly, this will sometimes place put points outside the model domain
%   (though I'm not yet sure why). The net result is that you have to
%   re-edit the grid SMS by deleting that particular node and recreating
%   the triangulation manually.
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2013-03-11 First version.
%   2013-03-19 Add optional minimum angle support.
%
%==========================================================================

subname = 'fix_inside_boundary';

global ftbverbose
if ftbverbose
    fprintf('\n'); fprintf(['begin : ' subname '\n']);
end

% Check the inputs
if length(x) ~= length(y)
    error('Size of inputs (x, y) do not match.')
else
    % Set the number of points in the boundary
    np = length(x(node_ids));
end

if nargin == 4
    minAng = true;
else
    minAng = false;
end

normx = nan(np, 1);
normy = nan(np, 1);

% Order the boundary points in clockwise order.
[boundx, boundy] = poly2cw(x(node_ids), y(node_ids));

% Create an array for the list of indices which have been moved which are
% adjacent to the open boundary.
seen = nan(np, 1);

% Create the output arrays from the input data. These values will be
% adjusted for the nodes adjacent to the open boundary defined by node_ids.
adjx = x;
adjy = y;

for pp = 1:np
    % For each node, find the two closest nodes to it and use those as the
    % adjacent vectors. Doing this is slower than just assuming the nodes
    % are provided in a sorted order, but means we don't have to worry too
    % much about any irregularities from the sorting above (poly2cw).
    [~, idx] = sort(sqrt((boundx(pp) - boundx).^2 + (boundy(pp) - boundy).^2));
    
    % Get the coordinates of the two nearest points (skip the first closest
    % because that's the current point).
    [px1, py1] = deal(boundx(idx(2)), boundy(idx(2)));
    [px2, py2] = deal(boundx(idx(3)), boundy(idx(3)));
    
    % Find the length of the edges of the triangle formed from the three
    % points we're currently considering. 1 and 2 are those adjacent and 3
    % is the remaining side.
    ln1 = sqrt((boundx(pp) - px1)^2 + (boundy(pp) - py1)^2);
    ln2 = sqrt((boundx(pp) - px2)^2 + (boundy(pp) - py2)^2);
    ln3 = sqrt((px1 - px2)^2 + (py1 - py2)^2);
    
    % Find the angle between the two element edges and the current node
    % (cosine rule). Use the real component only for cases where the three
    % points lie on a straight line (in which case the angle should be
    % 180 degrees).
    ang1 = real(acosd((ln1^2 + ln2^2 - ln3^2) / (2 * ln1 * ln2)));
    ang1b = ang1 / 2; % bisect the angle

    % Check if we've been given a threshold minimum angle and skip this
    % open boundary point if we have.
    if minAng
        if ang1 < ang_thresh
            continue
        end
    end

    % Find the angle to the horizontal for the current node and one of the
    % other points.
    ang2 = atan2((py1 - boundy(pp)), (px1 - boundx(pp))) * (180 / pi);
    
    % Find the difference between the two.
    ang3 = ang2 - ang1b;
    
    % Now get the mean length of the two closest element edges and use that
    % to create the new point inside the boundary. Scale it to 90% of the
    % value to make a cleaner transition when we're shifting points below.
    ml = 0.9 * mean([ln1, ln2]);
    dy = ml * sind(ang3);
    dx = ml * cosd(ang3);

    % Add the offsets to the current node to get the new node's position.
    [xx(1), yy(1)] = deal(boundx(pp) + dx, boundy(pp) + dy);
    [xx(2), yy(2)] = deal(boundx(pp) - dx, boundy(pp) - dy);

    % Check which of the two sets above is inside the polygon defined by
    % the open boundary.
    if inpolygon(xx(1), yy(1), boundx, boundy) == 1
        [normx(pp, 1), normy(pp, 1)] = deal(xx(1), yy(1));
    elseif inpolygon(xx(2), yy(2), boundx, boundy) == 1
        [normx(pp, 1), normy(pp, 1)] = deal(xx(2), yy(2));
    else
        warning('Both versions of the calculated point are outside the model domain. Skipping.')
        continue
    end
        
    % OK, so now we have a new point approximately orthogonal to the
    % current node, we can use its position to find the nearest existing
    % node inside the domain and replace its coordinates with the new
    % node's.
    [~, idx2] = sort(sqrt((normx(pp, 1) - x).^2 + (normy(pp, 1) - y).^2));
    
    % We need to check we haven't seen this node before (in 'seen') and
    % also that it's not an open boundary point.
    c = 1;
    while true
        if ismember(idx2(c), node_ids) || ismember(idx2(c), seen)
            % Keep going until we find one that's not been used before.
            c = c + 1;
        else
            break
        end
    end
    % Append to the list of coordinates we've seen so we don't move the
    % same node twice.
    seen(pp) = idx2(c);

    % Replace the coordinates.
    adjx(idx2(c)) = normx(pp, 1);
    adjy(idx2(c)) = normy(pp, 1);

end

if ftbverbose
    fprintf('end   : %s\n', subname)
end

% Do a figure to see the effect
% close all
% h = figure(1);
% 
% % Original node positions
% plot(x, y, 'o', 'MarkerEdgeColor', [0.9, 0.1, 0.1], ...
%     'MarkerFaceColor', [1, 0.1, 0.1])
% hold on
% axis('equal', 'tight')
% % Adjusted node positions
% plot(adjx, adjy, 'o', 'MarkerEdgeColor', [0.1, 0.1, 0.8], ...
%     'MarkerFaceColor', [0.1, 0.1, 0.9])
% 
% % Original triangulation
% patch('Vertices', [x, y], 'Faces', Mobj.tri, 'CData', [], ...
% 'edgecolor', [0.9, 0.1, 0.1], 'FaceColor', 'w', 'LineWidth', 1);
% % New triangulation
% patch('Vertices', [adjx, adjy], 'Faces', Mobj.tri, 'CData', [], ...
% 'edgecolor', [0.1, 0.1, 0.8], 'FaceColor', 'w', 'LineWidth', 1);
% 
% % Open boundary nodes
% % plot(x(node_ids), y(node_ids), 'ko', 'MarkerFaceColor', 'k')
% 
% legend('Original nodes', 'Adjusted nodes')
% legend('BoxOff')
% 
% xlabel('Eastings')
% ylabel('Northings')
