function Mobj = find_boundary_elements(Mobj)
% Find the elements which fall along the boundary.
% 
% Mobj = find_boundary_elements(Mobj)
% 
% DESCRIPTION:
%   Find the elements which are bounded by the open boundaries described by
%   the nodes in Mobj.read_obc_nodes.
% 
% INPUT:
%   Mobj - required fields:
%           - read_obc_nodes
%           - obc_nodes
%           - tri
% 
% OUTPUT:
%   Mobj - new field of a cell array read_obc_elements which contains the
%          IDs of the elements which fall on the model open boundaries and
%          nObcElements which is the total number of boundary elements
%          along each boundary.
% 
% NOTES:
%   This will be pretty slow if your unstructured grid has an enormous
%   number of elements in it (it loops through every element and compares
%   against the boundary nodes). I'm sure there's a quicker way, so feel
%   free to have at it.
% 
% EXAMPLE USAGE:
%   Mobj = find_boundary_elements(Mobj)
% 
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
% 
% Revision history:
%   2013-02-26 First version.
%   2013-02-28 Add new field to the output (total number of boundary
%   elements as nObcElements).
% 
%==========================================================================

subname = 'find_boundary_elements';

global ftbverbose
if ftbverbose
    fprintf('\n')
    fprintf('begin : %s\n', subname)
end

ne = length(Mobj.tri); % number of elements
nb = length(Mobj.read_obc_nodes); % number of boundaries

obc_elems = cell(nb, 1);
nObcElements = nan(nb, 1);

for i = 1:nb

    % Do the current boundary's nodes
    nodeIDs = Mobj.obc_nodes(i, Mobj.obc_nodes(i, :) ~= 0);

    f = 0;
    
    for ttt = 1:ne
        tri = Mobj.tri(ttt, :);
        C = intersect(tri, nodeIDs);
        % Only those with a face along the boundary count (i.e. two nodes
        % on the boundary), particularly for the mean flow.
        if numel(C) == 2
            f = f + 1; % increment the found counter
            obc_elems{i}(f) = ttt;
            if ftbverbose
                fprintf('Found boundary element ID %i\n', ttt)
            end
        end
    end
    nObcElements(i) = numel(obc_elems{i}(:));
end

Mobj.read_obc_elements = obc_elems;
Mobj.nObcElements = nObcElements;

% Check it's worked for the first model boundary.
% xc = nodes2elems(Mobj.x, Mobj);
% yc = nodes2elems(Mobj.y, Mobj);
% figure(1)
% clf
% plot(Mobj.x, Mobj.y, 'r.', xc, yc, 'ko')
% hold on
% plot(xc(Mobj.read_obc_elements{1}), yc(Mobj.read_obc_elements{1}), 'gx')
% axis('equal', 'tight')

if ftbverbose
    fprintf('end   : %s \n', subname)
end
