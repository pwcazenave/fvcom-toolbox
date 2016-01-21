function [Nested]=find_nesting_region(conf,Mobj)
% Creates a nesting structure array for direct/indirect or weighted nesting
%
% function [Nested]=find_nesting_region(conf,M)
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
%              make_input.m type of script. Minimum fields are
%              Nesting type (1, 2 == direct nesting, 3 == weighted)
%              conf.type = 3;
%              If we're doing type 3, we can specify the number of levels of nested
%              boundaries to use. The minimum valid value is 1. For Indirect or direct
%              nesting use 1
%              conf.levels = 5;
%              conf.power = determines drop of weights from 1 [0 is linear, 1-4 is 1/conf.levels.^conf.power]
%  M        = Mesh object
%
% OUTPUT:
%  Nested  Mesh object with added nesting variables .
%
% EXAMPLE USAGE:
%   conf.type = type of nesting  [1, 2 == direct nesting, 3 == weighted]
%   conf.levels = number of boundary bands to use for weighted option
%   conf.power = determines drop of weights from 1 [0 is linear, 1-4 is 1/conf.levels.^conf.power]
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
%       in setup_metrics.m
%
%   the global variable ftbverbose shows information at run time as well as
%   generating a figure with the location of the nesting nodes and elements
%
%   [Nested]=find_nesting_region(conf,Mobj)
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
%   Torres by Pierre
%   script.
%   2016-01-19 Updated to a stand alone function and general tidy up
%   general tidy up.
%
%==========================================================================


dump = dbstack;
subname = dump.name;
clear dump

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

M=Mobj;

TR = triangulation(M.tri, [M.x, M.y]);
et = cell(conf.levels, 1);
M.weight_cell = cell(1,length(M.nObcNodes)*conf.levels);
for oo=1:length(M.read_obc_nodes)
    M.weight_node{oo} = ones(1, M.nObcNodes(oo));
end
nn=0;
if ftbverbose
    figure,cla
    patch('Vertices', [M.x, M.y], 'Faces', M.tri, ...
        'Cdata', -M.h, 'edgecolor', 'b', 'facecolor', 'interp');
end
%% The initial OB are first in the order. The nested OB
% follow

etdx=0;
for oo=1:length(Mobj.read_obc_nodes)
    for n = 1:conf.levels
        nn=nn+1;etdx=etdx+1;
        if (oo>1 && n == 1)
            M.read_obc_nodes{nn+1} = Mobj.read_obc_nodes{oo};
            M.obc_type(nn+1)=Mobj.obc_type(oo);
            M.nObcNodes(nn+1) = length(Mobj.read_obc_nodes{oo});
            M.obc_nodes(nn+1, 1:Mobj.nObcNodes(oo)) = Mobj.read_obc_nodes{oo};
            
            ti = vertexAttachments(TR, double(M.read_obc_nodes{nn+1})');
            et{etdx} = setdiff(unique([ti{:}]),[et{1:end}]);
            nn=nn+1;
            
            M.read_obc_nodes{nn+1} = int32(setdiff(unique(M.tri(et{etdx}, :)), ...
                [M.read_obc_nodes{1:nn}]))';
            M.obc_type(nn+1)=M.obc_type(nn);
            M.nObs =  M.nObs+1;
            M.nObcNodes(nn+1) = length(M.read_obc_nodes{nn+1});
            M.obc_nodes(nn+1, 1:M.nObcNodes(nn+1)) = M.read_obc_nodes{nn+1};
        else
            
            ti = vertexAttachments(TR, double(M.read_obc_nodes{nn})');
            et{etdx} = setdiff(unique([ti{:}]),[et{1:end}]);
            M.read_obc_nodes{nn+1} = int32(setdiff(unique(M.tri(et{etdx}, :)), ...
                [M.read_obc_nodes{1:nn}]))';
            M.obc_type(nn+1)=M.obc_type(nn);
            M.nObs =  M.nObs+1;
            M.nObcNodes(nn+1) = length(M.read_obc_nodes{nn+1});
            M.obc_nodes(nn+1, 1:M.nObcNodes(nn+1)) = M.read_obc_nodes{nn+1};
            
        end
        % Plot nesting region.
        if ftbverbose
            if n == 1 | nn==1
                axis('equal', 'tight')
                colormap('gray')
                hold on
                plot(M.x(M.read_obc_nodes{nn}), M.y(M.read_obc_nodes{nn}), 'wo')
                plot(M.x(M.read_obc_nodes{nn+1}), M.y(M.read_obc_nodes{nn+1}), 'ro')
                plot(M.xc(et{etdx}), M.yc(et{etdx}), 'wx')
            else
                plot(M.x(M.read_obc_nodes{nn+1}), M.y(M.read_obc_nodes{nn+1}), 'ro')
                plot(M.xc(et{etdx}), M.yc(et{etdx}), 'wx')
            end
        end
        
    end
end
if conf.type == 3
    % probably best to have them as a linear decrease?
    if conf.power ==0
        % one should be the outermost boundary level which is why the vectors
        % are flipped
        weights_nodes = fliplr((1:conf.levels)./conf.levels);
        weights_nodes(end+1) =0;
        
        weights_elems = fliplr((1:conf.levels-1)./conf.levels-1);
        weights_elems(end+1) =0;
        
    else
        weights_nodes = 1:conf.levels+1;
        weights_nodes = 1./weights_nodes.^conf.power;
        weights_elems = 1:conf.levels;
        weights_elems = 1./weights_elems.^conf.power;
    end
    
    nn=0;etdx=0;
    for oo=1:length(Mobj.read_obc_nodes)
        nn=nn+1;
        
        M.weight_node{nn} = repmat(weights_nodes(1), ...
            1, M.nObcNodes(nn));
        for n = 1:conf.levels
            etdx=etdx+1;
            M.weight_node{nn+1} = repmat(weights_nodes(n+1), ...
                1, M.nObcNodes(nn+1));
            M.weight_cell{etdx} = repmat(weights_elems(n), ...
                1, length(et{etdx}), 1);
            nn=nn+1;
        end
    end
end
Nested=M;
Nested.read_obc_elems=et;
return

