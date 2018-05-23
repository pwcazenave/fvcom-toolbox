function Nested = find_nesting_region(conf, Mobj)
% Creates a nesting structure array for direct/indirect or weighted nesting
%
% function Nested = find_nesting_region(conf, Mobj)
%
% DESCRIPTION:
%   Uses the Mesh object Mobj and the conf variable to search for the nodes
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
%   Nested.obc_nodes = matrix with node indices. Each row is a boundary level.
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
%   2016-12-14 Updated the help. Also disabled the plot to ease automated
%   runs.
%   2016-12-22 Fairly major rewrite to make things clearer and less prone
%   to subtle bugs.
%   2017-02-16 Fix for direct nesting (no weights needed).
%   2017-02-20 Fix non-critical bug which added empty cells for the nest
%   elements. Also make the node and element weight values match one another.
%
%==========================================================================

[~, subname] = fileparts(mfilename('fullpath'));

global ftbverbose
if isempty(ftbverbose)
    ftbverbose=false;
end
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

TR = triangulation(Mobj.tri, [Mobj.x, Mobj.y]);

Nested = Mobj;
Nested.nObs = 0; % number of nodal levels is incremented for each level.

% Make cell arrays to store the element IDs for each nested level as well
% as for the weights on the nodes and elements.
Nested.read_obc_nodes = cell(0);
Nested.read_obc_elems = cell(0);
if any(conf.Nested_type == 3)
    Nested.weight_cell = cell(0);
    Nested.weight_node = cell(0);
end

% if ftbverbose
%     figure(1)
%     clf
%     triplot(Nested.tri, Nested.x, Nested.y)
%     axis('equal', 'tight')
%     hold on
% end

% Indices for the output cell arrays which are incremented for each nest
% level and each open boundary.
cumulative_node_idx = 1;
cumulative_elem_idx = 1;
for obc_idx = 1:Mobj.nObs
    % Generate the weights for the elements and nodes.
    if conf.Nested_type(obc_idx) == 3
        if conf.power == 0
            weights_nodes = fliplr((1:conf.levels(obc_idx) + 1)./(conf.levels(obc_idx) + 1));
        else
            weights_nodes = 1:conf.levels(obc_idx) + 1;
            weights_nodes = 1./weights_nodes.^conf.power;
        end
        % Use the same weights as the nodes, but drop the last level
        % (elements is always 1 smaller).
        weights_elems = weights_nodes(1:end - 1);
    end

    % Save the original open boundary nodes into the nested struct (Nested).
    Nested.read_obc_nodes{cumulative_node_idx} = Mobj.read_obc_nodes{obc_idx};

    % Given the current open boundary, find the elements connected to it
    % and give them some weights.
    ti = vertexAttachments(TR, double(Mobj.read_obc_nodes{obc_idx})');
    Nested.read_obc_elems{cumulative_elem_idx} = unique([ti{:}]);

    % Save the weights into the nested struct (Nested).
    if conf.Nested_type(obc_idx) ~= 1
        Nested.weight_node{cumulative_node_idx} = repmat(weights_nodes(1), 1, length(Nested.read_obc_nodes{cumulative_node_idx}));
        Nested.weight_cell{cumulative_elem_idx} = repmat(weights_elems(1), 1, length(Nested.read_obc_elems{cumulative_elem_idx}));
    end

    % Also save the type of open boundary we've got and update the open
    % boundary counter and number of open boundary nodes.
    Nested.nObcNodes(cumulative_node_idx) = length(Nested.read_obc_nodes{cumulative_node_idx});
    Nested.nObs = Nested.nObs + 1;
    Nested.obc_type(cumulative_node_idx) = conf.Nested_type(obc_idx);

    if ftbverbose && conf.Nested_type(obc_idx) ~= 1
%         figure(1)
%         scatter(Nested.x(Nested.read_obc_nodes{cumulative_node_idx}), Nested.y(Nested.read_obc_nodes{cumulative_node_idx}), 20, Nested.weight_node{cumulative_node_idx}, 'filled')
%         scatter(Nested.xc(Nested.read_obc_elems{cumulative_elem_idx}), Nested.yc(Nested.read_obc_elems{cumulative_elem_idx}), 20, Nested.weight_cell{cumulative_elem_idx}, 'filled')
        fprintf('Original open boundary %d\n', obc_idx)
    end

    % Now we have the original open boundary and the elements connected to
    % it we can move through the levels specified in conf.levels(obc_idx)
    % and repeat the process. Bump the cumulative counters accordingly.
    cumulative_node_idx = cumulative_node_idx + 1;
    cumulative_elem_idx = cumulative_elem_idx + 1;

    for lev = 1:conf.levels(obc_idx)
        % Find the nodes and elements for this level and assign their
        % weights. Use the most recent data in Nested.read_obc_nodes as the
        % anchor from which to work.
        Nested.read_obc_nodes{cumulative_node_idx} = int32(setdiff(unique(Nested.tri(Nested.read_obc_elems{cumulative_elem_idx - 1}, :)), ...
            [Nested.read_obc_nodes{1:cumulative_node_idx - 1}]))';
        ti = vertexAttachments(TR, double(Nested.read_obc_nodes{cumulative_node_idx})');
        Nested.nObs = Nested.nObs + 1;
        Nested.obc_type(cumulative_node_idx) = conf.Nested_type(obc_idx);
        Nested.nObcNodes(cumulative_node_idx) = length(Nested.read_obc_nodes{cumulative_node_idx});

        if conf.Nested_type(obc_idx) ~= 1
            Nested.weight_node{cumulative_node_idx} = repmat(weights_nodes(lev + 1), 1, length(Nested.read_obc_nodes{cumulative_node_idx}));
        end
        if lev ~= conf.levels(obc_idx)
            Nested.read_obc_elems{cumulative_elem_idx} = setdiff(unique([ti{:}]), [Nested.read_obc_elems{:}]);
            if conf.Nested_type(obc_idx) ~= 1
                Nested.weight_cell{cumulative_elem_idx} = repmat(weights_elems(lev + 1), 1, length(Nested.read_obc_elems{cumulative_elem_idx}));
            end
        end

        if ftbverbose && conf.Nested_type(obc_idx) ~= 1
%             figure(1)
%             scatter(Nested.x(Nested.read_obc_nodes{cumulative_node_idx}), Nested.y(Nested.read_obc_nodes{cumulative_node_idx}), 20, Nested.weight_node{cumulative_node_idx}, 'filled')
%             if lev ~= conf.levels(obc_idx)
%                 scatter(Nested.xc(Nested.read_obc_elems{cumulative_elem_idx}), Nested.yc(Nested.read_obc_elems{cumulative_elem_idx}), 20, Nested.weight_cell{cumulative_elem_idx}, 'filled')
%             end
            fprintf('Nested level %d\n', lev)
        end

        % Bump the node and element cumulative counters so the next loop
        % dumps everything into the right position in the cell arrays.
        cumulative_node_idx = cumulative_node_idx + 1;
        if lev ~= conf.levels(obc_idx)
            cumulative_elem_idx = cumulative_elem_idx + 1;
        end
    end
    % Check if all possible elements have been correctly identified
    % colapse all nodes into an array and extract all elements connected to
    % those nodes
    % find the elements that are attached to 3 nodes in the nesting region
    %      
    candidate= find(all(ismember(Mobj.tri,[Nested.read_obc_nodes{end}]),2));
    % test if it is different from existing elements
    [truecandidate]=setdiff(candidate,[Nested.read_obc_elems{:}]);
    % add to existing list. Catenate will ignore an empty result from
    % setdiff
    if ~isempty(truecandidate)
        for dd=1:length(truecandidate)
    Nested.read_obc_elems{end} = cat(2,Nested.read_obc_elems{end},truecandidate(dd));
    Nested.weight_cell{end} = cat(2,Nested.weight_cell{end},Nested.weight_cell{end}(end));
    end
    end
    if ftbverbose
        fprintf('\n')
    end
end

% Update the clunky obc_nodes array with the new node IDs from
% Nested.read_obc_nodes.
for nidx = 1:length(Nested.read_obc_nodes)
    Nested.obc_nodes(nidx, 1:length(Nested.read_obc_nodes{nidx})) = Nested.read_obc_nodes{nidx};
end

if ftbverbose
    % figure(1)
    % colorbar
    % title('Nest weights')
    fprintf('end   : %s \n', subname)
end
