function Nested = find_nesting_region(conf, Mobj)
% Creates a nesting structure array for direct/indirect or weighted nesting
%
% function Nested = find_nesting_region(conf,M)
%
% DESCRIPTION:
%   Uses the Mesh object M and the conf variable to search for the nodes
%   and elements originating from the open boundaries into de domain
%   interior for a specified number of levels.
%
% Optionally specify nesting type:
%   1/2: DIRECT/INDIRECT nesting:
%       - Full variables/no surface elevation respectively.
%   3:   RELAXATION nesting:
%       - Nesting with a relaxation method.
%
% INPUT:
%  conf     = struct whose field names are the model configuration options
%              as set by the user. These are generally assigned at the top of a
%              make_input.m type of script. Minimum fields are:
%              - Nested_type: Array of nesting types, one for each open boundary.
%              1 or 2 are direct nesting, 3 is weighted. If we're doing
%              type 3, we can specify the number of levels of nested
%              boundaries to use. The minimum valid value is 1. For
%              Indirect or direct nesting use 1
%              - levels: number of levels of nodes over which to relax the
%              boundary (if Nested_type is 3).
%              - power = determines drop of weights from 1. 0 is linear,
%              1-4 is 1/conf.levels.^conf.power.
%  Mobj     = Mesh object with the following fields:
%              - tri: Triangulation table as pupulated by read_sms_grid
%              - x: Node x coordinates (cartesian)
%              - y: Node y coordinates (cartesian)
%              - xc: element x coordinates (cartesian)
%              - yc: Element y coordinates (cartesian)
%              - nObcNodes: number of nodes as defined by SMS and read
%              through read_sms_grid
%              - read_obc_nodes = nodes indices at boundaries as read in
%              read_sms_grid.
%              - obc_type = Type of OBC as defined in mod_obcs.F [ values
%              of 1-10] and parsed to Mobj from conf by add_obc_nodes_list
%              - obc_nodes: matrix with node indices. Each row is a given
%              open boundary.
%              - nObs: number of open boundaries.
%
% OUTPUT:
%  Nested   = Mesh object with all the same fields as Mobj, plus the
%             following added and modified fields:
%              - read_obc_nodes: new nested boundary node IDs
%              - read_obc_elems: new nested boundary element IDs
%              - nObs: number of open boundaries (number of levels * number
%              - nObcNodes: number of nodes in each nested level
%              of original open boundaries)
%              - obc_type: the type for each nested boundary
%              - obc_nodes: the array-based nested boundary node IDs (for
%              backwards compatibility)
%              - weight_node: weights for each nested node boundary level
%              - weight_cell: weights for each nested element boundary level
%
% EXAMPLE USAGE:
%   conf.Nested_type = type of nesting [1, 2 == direct nesting, 3 == weighted]
%   conf.levels = number of boundary bands to use for weighted option
%   conf.power = determines drop of weights from 1 [0 is linear, anything
%   else is 1/conf.levels.^conf.power]
%   Mobj.tri = Triangulation table as pupulated by read_sms_grid
%   Mobj.x = Nodes x coordinate [we should make it possible to use lon instead]
%   Mobj.y = Nodes y coordinate [we should make it possible to use lat instead]
%   Mobj.nObcNodes = number of nodes as defined by SMS and read through read_sms_grid
%   Mobj.read_obc_nodes = nodes indices at boundaries as read in read_sms_grid
%   Mobj.obc_type = Type of OBC as defined in mod_obcs.F [ values of 1-10]
%       and parsed to Mobj from conf by add_obc_nodes_list
%   Mobj.xc = Nodes x coordinate [we should make it possible to use lonc instead]
%   Mobj.yc = Elements y coordinate [we should make it possible to use latc instead]
%   Mobj.obc_nodes = matrix with node indices. Each row is a boundary  level.
%   Mobj.nObs = total number of open boundary levels (I think this is set
%       in setup_metrics.m).
%
%   the global variable ftbverbose shows information at run time as well as
%   generating a figure with the location of the nesting nodes and elements
%
%   Nested = find_nesting_region(conf,Mobj)
%
%   Nested.nObcNodes = number of nodes in new nesting zone
%   Nested.read_obc_nodes = nodes indices in nesting zone
%   Nested.obc_type = Type of OBC as defined in mod_obcs.F [ values of 1-10]
%   Nested.obc_nodes = matrix with node indices. Each row is a boundary  level.
%   Nested.nObs = total number of open boundary levels
%   Nested.weight_node = weights for nodes if using weighted type [0-1] see
%       FVCOM manual for further info (Chapter 6 section 4 in version 3.2)
%   Nested.weight_cell = weights for elements if using weighted type [0-1]

% Author(s):
%   Ricardo Torres (Plymouth Marine Laboratory)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%   Darren Price (CH2MHill)
%   Hakeem Johnson (CH2MHill)
%
% Revision history:
%   2015-11-01 First version based on Hakeem and Darren code as provided to
%   Torres by Pierre.
%   2016-01-19 Updated to a stand alone function and general tidy up.
%   2016-12-14 Updated the help.
%
%==========================================================================


dump = dbstack;
subname = dump.name;
clear dump

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

M = Mobj;
M.nObcNodes = M.nObcNodes.*(~isnan(conf.levels./conf.levels));
nBC = sum(~isnan(conf.levels./conf.levels));
TR = triangulation(M.tri, [M.x, M.y]);
et = cell(sum(conf.levels), 1);
M.weight_cell = cell(1, nBC*sum(conf.levels));
for oo = 1:length(M.read_obc_nodes)
    M.weight_node{oo} = ones(1, M.nObcNodes(oo));
end
nn = 0;
if ftbverbose
    figure,cla
    patch('Vertices', [M.x, M.y], 'Faces', M.tri, ...
        'Cdata', -M.h, 'edgecolor', 'b', 'facecolor', 'interp');
end

% The initial boundary is first in the order. The nested boundaries follow.
etdx = 0;
for oo=1:length(Mobj.read_obc_nodes)
    if conf.Nested_type(oo) == 3
        % probably best to have them as a linear decrease?
        if conf.power ==0
            % one should be the outermost boundary level which is why the vectors
            % are flipped
            weights_nodes = fliplr((1:conf.levels(oo))./conf.levels(oo));
            weights_nodes(end+1) = 0;
            
            weights_elems = fliplr((1:conf.levels(oo)-1)./conf.levels(oo)-1);
            weights_elems(end+1) = 0;
            
        else
            weights_nodes = 1:conf.levels(oo)+1;
            weights_nodes = 1./weights_nodes.^conf.power;
            weights_elems = 1:conf.levels(oo);
            weights_elems = 1./weights_elems.^conf.power;
        end
        M.weight_node{oo} = repmat(weights_nodes(1), ...
            1, M.nObcNodes(oo));
        
    end
    for n = 1:conf.levels(oo)
        nn = nn + 1;
        etdx = etdx + 1;
        if (oo > 1 && n == 1)
            M.read_obc_nodes{nn + 1} = Mobj.read_obc_nodes{oo};
            M.obc_type(nn + 1)=conf.Nested_type(oo);
            M.nObcNodes(nn + 1) = length(Mobj.read_obc_nodes{oo});
            M.obc_nodes(nn + 1, 1:Mobj.nObcNodes(oo)) = Mobj.read_obc_nodes{oo};
            
            ti = vertexAttachments(TR, double(M.read_obc_nodes{nn + 1})');
            et{etdx} = setdiff(unique([ti{:}]),[et{1:end}]);
            nn=nn + 1;
            
            M.read_obc_nodes{nn + 1} = int32(setdiff(unique(M.tri(et{etdx}, :)), ...
                [M.read_obc_nodes{1:nn}]))';
            M.obc_type(nn + 1) = conf.Nested_type(oo);
            M.nObs =  M.nObs+1;
            M.nObcNodes(nn + 1) = length(M.read_obc_nodes{nn + 1});
            M.obc_nodes(nn + 1, 1:M.nObcNodes(nn + 1)) = M.read_obc_nodes{nn + 1};
            
        else
            
            ti = vertexAttachments(TR, double(M.read_obc_nodes{nn})');
            et{etdx} = setdiff(unique([ti{:}]),[et{1:end}]);
            M.read_obc_nodes{nn + 1} = int32(setdiff(unique(M.tri(et{etdx}, :)), ...
                [M.read_obc_nodes{1:nn}]))';
            M.obc_type(nn + 1)=M.obc_type(nn);
            M.nObs =  M.nObs+1;
            M.nObcNodes(nn + 1) = length(M.read_obc_nodes{nn + 1});
            M.obc_nodes(nn + 1, 1:M.nObcNodes(nn + 1)) = M.read_obc_nodes{nn + 1};
            
        end
        % Plot nesting region.
        if ftbverbose
            if n == 1 || nn == 1
                axis('equal', 'tight')
                colormap('gray')
                hold on
                plot(M.x(M.read_obc_nodes{nn}), M.y(M.read_obc_nodes{nn}), 'wo')
                plot(M.x(M.read_obc_nodes{nn + 1}), M.y(M.read_obc_nodes{nn + 1}), 'ro')
                plot(M.xc(et{etdx}), M.yc(et{etdx}), 'wx')
            else
                plot(M.x(M.read_obc_nodes{nn + 1}), M.y(M.read_obc_nodes{nn + 1}), 'ro')
                plot(M.xc(et{etdx}), M.yc(et{etdx}), 'wx')
            end
        end
        
        
        if conf.Nested_type(oo) == 3
            M.weight_node{nn + 1} = repmat(weights_nodes(n+1), ...
                1, M.nObcNodes(nn + 1));
            M.weight_cell{etdx} = repmat(weights_elems(n), ...
                1, length(et{etdx}), 1);
        end
    end
end

Nested = M;
Nested.read_obc_elems = et;

if ftbverbose
    fprintf('end   : %s \n', subname)
end

return
