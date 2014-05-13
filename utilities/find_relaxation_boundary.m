function Mobj = find_relaxation_boundary(Mobj)
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

nb = length(Mobj.read_obc_nodes); % number of boundaries
bc_width=Mobj.relax_bc_Nnodes;
obc_elems = cell(nb, 1);
nObcElements = nan(nb, 1);


for i = 1:nb
    
    % Do the current boundary's nodes
    nodeIDs = Mobj.obc_nodes(i, Mobj.obc_nodes(i, :) ~= 0);
    [C1,~] = ismember(Mobj.tri(:,1),nodeIDs(:), 'rows');
    [C2,~] = ismember(Mobj.tri(:,2), nodeIDs(:),'rows');
    [C3,~] = ismember(Mobj.tri(:,3), nodeIDs(:),'rows');
    obc_elems{i}= unique([find(C1);find(C2);find(C3)]);
    nObcElements(i) = numel(obc_elems{i}(:));
    
end
Mobj.relaxBC_nodes={reshape([Mobj.read_obc_nodes{:}],[],1)};
Mobj.relaxBC_elems={reshape([obc_elems{:}],[],1)};

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


% nodeIDs = Mobj.tri(Mobj.relaxBC_elems{bb},:);
% nodeIDs=unique(nodeIDs(:));
% C1 = setdiff(nodeIDs(:),...
%     cat(1,Mobj.relaxBC_nodes{1:bb}), 'rows');
% Mobj.relaxBC_nodes(bb+1)={C1};
Mobj.relaxnBC_elems=length(cat(1,Mobj.relaxBC_elems{:}));
Mobj.relaxnBC_nodes=length(cat(1,Mobj.relaxBC_nodes{:}));
% % Check it's worked for the first model boundary.
% xc = nodes2elems(Mobj.x, Mobj);
% yc = nodes2elems(Mobj.y, Mobj);
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
