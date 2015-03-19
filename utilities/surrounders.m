function nodes = surrounders(n, triangles)
% Return the IDs of the nodes surrounding node number `n'.
%
% INPUTS:
%   n
%       Node ID around which to find the connected nodes.
%   triangles
%       Triangulation matrix to find the connected nodes.
%
% OUTPUTS:
%   surroundingidx
%       Indices of the surrounding nodes.
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2015-02-19 First version based on the PyFVCOM.grid_tools surrounders
%   function.
%
% Notes:
%
% Check it works with:
% [x, y] = meshgrid(1:25, 100:125);
% x = x(:) + randn(numel(x), 1) * 0.1;
% y = y(:) + randn(numel(y), 1) * 0.1;
% tri = delaunay([x, y]);
% for n = linspace(1, length(x) - 1, 5)
% 	  aa = surrounders(n, tri)
% 	  figure
% 	  triplot(x, y, tri)
% 	  plot(x(n), y(n), 'ro')
% 	  plot(x(aa), y(aa), 'ko')
% 	  xlim(min(x(aa)) - 1, max(x(aa)) + 1)
% 	  ylim(min(y(aa)) - 1, max(y(aa)) + 1)
% 	  legend
% end

eidx = max((abs(triangles - n) == 0), [], 2);
nodes = unique(triangles(triangles(eidx) ~= n, :));
