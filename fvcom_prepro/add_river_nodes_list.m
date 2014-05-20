function [Mobj] = add_river_nodes_list(Mobj,Nlist,RiverName,plotFig)

% Add a set of river nodes comprising a single river to Mesh structure  
% Using a set of user-defined nodes
%
% [Mobj] = add_river_nodes(Mobj,Nlist,RiverName)
%
% DESCRIPTION:
%    Select using ginput the set of nodes comprising a river
%
% INPUT
%    Mobj = Matlab mesh object
%    RiverName = Name of the River
%    plotFig = [optional] show a figure of the mesh (1 = yes)
%
% OUTPUT:
%    Mobj = Matlab mesh object with an additional river nodelist
%
% EXAMPLE USAGE
%    Mobj = add_river_nodes(Mobj, [146, 3004], 'Potomac')
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Karen Amoudry (National Oceanography Centre, Liverpool)
%
% Note:
%    Uses ginput2 which allows zooming before selecting points and displays
%    clicked points realtime
%
% Revision history
%    2013-01-02 KJA bug fix: amended usage of 'unique' to prevent it from
%    sorting the values it returns. Amended by Pierre to support pre-2012
%    versions of MATLAB whilst giving the same result.
%    2013-10-22 KJA: added capability to turn off figures (copied from
%    Pierre's adaptation to add_obc_nodes_list.m)
%    2014-05-20 Set boolean flag to true to indicate rivers.
%   
%==========================================================================
subname = 'add_river_nodes_list';
global ftbverbose
if ftbverbose
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end

% Do we want a figure showing how we're getting along?
if nargin == 3
    plotFig = 0;
end

%--------------------------------------------------------------------------
% Get a unique list and make sure they are in the range of node numbers 
%--------------------------------------------------------------------------
% Make this works in versions of MATLAB older than 2012a (newer versions
% can just use unique(A, 'stable'), but checking versions is a pain).
[~, Nidx] = unique(Nlist);
Nlist = Nlist(sort(Nidx));

if max(Nlist) > Mobj.nVerts
    fprintf('your river node number(s) exceed the total number of nodes in the domain\n');
    fprintf('stop screwing around\n');
    error('stopping...\n')
end

%--------------------------------------------------------------------------
% Plot the mesh 
%--------------------------------------------------------------------------
if plotFig == 1
    if strcmpi(Mobj.nativeCoords(1:3), 'car')
        x = Mobj.x;
        y = Mobj.y;
    else
        x = Mobj.lon;
        y = Mobj.lat;
    end
    
    figure
    patch('Vertices', [x,y], 'Faces', Mobj.tri,...
        'Cdata', Mobj.h, 'edgecolor', 'k', 'facecolor', 'interp');
    hold on
    
    plot(x(Nlist), y(Nlist), 'ro')
    title('river nodes')
end

% add to mesh object
npts = numel(Nlist);
Mobj.nRivers = Mobj.nRivers + 1;
Mobj.nRivNodes(Mobj.nRivers) = npts;
Mobj.riv_nodes(Mobj.nRivers, 1:npts) = Nlist;
Mobj.riv_name{Mobj.nRivers} = RiverName;

Mobj.have_rivers = true;

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end

