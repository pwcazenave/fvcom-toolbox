function Mobj = find_relaxation_boundary(Mobj)
% Find the elements which fall along the boundary for a nested
% configuration.
%
% Mobj = find_boundary_elements(Mobj)
%
% DESCRIPTION:
%   Find the elements which are bounded by the open boundaries described by
%   the nodes in Mobj.read_obc_nodes and which go to a depth of
%   Mobj.relax_bc_Nnodes.
%
% INPUT:
%   Mobj - required fields:
%           - read_obc_nodes - cell array of open boundary node IDs.
%           - tri - mesh triangulation
%           - relax_bc_Nnodes - array (length is Mobj.nObs) of number of
%           elements from the open boundary over which to relax the nested
%           inputs.
%
% OUTPUT:
%   Mobj - mesh object with the following new fields:
%           - relaxBC_nodes = node IDs for the relaxed region
%           - relaxBC_elems = element IDs for the relaxed region
%           - relaxnBC_nodes = number of nodes in the relaxed region
%           - relaxnBC_elems = number of elements in the relaxed region
%
% EXAMPLE USAGE:
%   Mobj = find_relaxation_boundary(Mobj)
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%   Ricardo Torres (Plymouth Marine Laboratory)
%
% Revision history:
%   2016-12-15 Add support for varying depth of nested region over each
%   open boundary. Also update help to actually refer to what this function
%   does.
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%   Ricardo Torres (Plymouth Marine Laboratory)
%
% Revision history:
%   2016-12-15 Update help to actually refer to what this function does.
%
%==========================================================================

subname = 'find_boundary_elements';

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

nb = length(Mobj.read_obc_nodes); % number of boundaries
bc_width = Mobj.relax_bc_Nnodes;
obc_elems = cell(nb, 1);
nObcElements = nan(nb, 1);

for i = 1:nb
    
    % Do the current boundary's nodes
    nodeIDs = Mobj.read_obc_nodes{i};
    [C1,~] = ismember(Mobj.tri(:,1), nodeIDs(:), 'rows');
    [C2,~] = ismember(Mobj.tri(:,2), nodeIDs(:), 'rows');
    [C3,~] = ismember(Mobj.tri(:,3), nodeIDs(:), 'rows');
    obc_elems{i} = unique([find(C1); find(C2); find(C3)]);
    if iscolumn(obc_elems{i})
        obc_elems{i} = obc_elems{i}';
    end
    nObcElements(i) = numel(obc_elems{i}(:));
    
end
Mobj.relaxBC_nodes = {[Mobj.read_obc_nodes{:}]};
Mobj.relaxBC_elems = {[obc_elems{:}]};

for bb=2:bc_width
    nodeIDs = Mobj.tri(Mobj.relaxBC_elems{bb-1},:);
    nodeIDs=unique(nodeIDs(:));
    C1 = setdiff(nodeIDs(:),...
        cat(1,Mobj.relaxBC_nodes{1:bb-1}), 'rows');
    Mobj.relaxBC_nodes(bb)={C1};
    [C1,~] = ismember(Mobj.tri(:,1),nodeIDs(:), 'rows');
    [C2,~] = ismember(Mobj.tri(:,2), nodeIDs(:),'rows');
    [C3,~] = ismember(Mobj.tri(:,3), nodeIDs(:),'rows');
    C1 = setdiff(unique([find(C1);find(C2);find(C3)]),...
        cat(1,Mobj.relaxBC_elems{1:bb-1}), 'rows');
    Mobj.relaxBC_elems(bb)={C1};
end


% figure(1)
% clf
%     triplot(Mobj.tri,Mobj.x,Mobj.y,'k');
% hold on
% 
% plot(Mobj.x( Mobj.relaxBC_nodes{1}), Mobj.y( Mobj.relaxBC_nodes{1}), 'r.')
% plot(Mobj.x( Mobj.relaxBC_nodes{2}), Mobj.y( Mobj.relaxBC_nodes{2}), 'r.')
% plot(Mobj.x( Mobj.relaxBC_nodes{3}), Mobj.y( Mobj.relaxBC_nodes{3}), 'kx')
% plot(Mobj.x( Mobj.relaxBC_nodes{4}), Mobj.y( Mobj.relaxBC_nodes{4}), 'rx')
% 
% plot(xc(Mobj.relaxBC_elems{3}), yc(Mobj.relaxBC_elems{3}), 'kx')
% axis('equal', 'tight')

if ftbverbose
    fprintf('end   : %s \n', subname)
end
