function Mobj = add_obc_nodes_list(Mobj, Nodelist, ObcName, ObcType, plotFig)
% Add a set of obc nodes comprising a single obc boundary to Mesh structure
% using a list of nodes.
%
% Mobj = add_obc_nodes_list(Mobj, Nodelist, ObcName, ObcType)
%
% DESCRIPTION:
%    Add a set of open boundary nodes for a given boundary using a list of
%    nodes.
%
% INPUT
%    Mobj = Matlab mesh object with fields:
%       - nVerts - number of nodes in the domain
%       - nativeCoords - coordinates in which the model runs (only for
%       plotting the figure).
%       - x, y, lon, lat - mesh node coordinates (either cartesian or
%       spherical) (only for plotting the figure).
%       - tri - model grid triangulation (only for plotting the figure)
%       - h - model grid depths (only for plotting the figure)
%    Nodelist = List of nodes
%    ObcName = Name of the Open Boundary
%    ObcType = FVCOM Flag for OBC Type
%    plotFig = [optional] show a figure of the mesh (1 = yes)
%
% OUTPUT:
%    Mobj = Matlab mesh object with an additional obc nodelist
%
% EXAMPLE USAGE
%    Nodelist = 1:100;
%    Mobj = add_obc_nodes_list(Mobj, Nodelist, 'OpenOcean')
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Karen Amoudry (National Oceanography Centre, Liverpool)
%
% Revision history:
%    2012-11-26 Add ability to turn off the figures.
%    2013-01-02 KJA bug fix: amended usage of 'unique' in line 53 to
%    prevent it from sorting the values it returns. Amended by Pierre to
%    support pre-2012 versions of MATLAB whilst giving the same result.
%    2015-02-23 Output number of nodes if the verbose flag is set.
%    2017-08-31 Update the help to clarify what's needed.
%
%==========================================================================
[~, subname] = fileparts(mfilename('fullpath'));
global ftbverbose
if ftbverbose
  fprintf('\nbegin : %s\n', subname)
end

% Do we want a figure showing how we're getting along?
if nargin == 4
    plotFig = 0;
end

%--------------------------------------------------------------------------
% Get a unique list and make sure they are in the range of node numbers
%--------------------------------------------------------------------------
% Make this works in versions of MATLAB older than 2012a (newer versions
% can just use unique(A, 'stable'), but checking versions is a pain).
[~, Nidx] = unique(Nodelist);
Nodelist = Nodelist(sort(Nidx));

assert(max(Nodelist) <= Mobj.nVerts, 'Your open boundary node number (%d) exceeds the total number of nodes in the domain (%d)', max(Nodelist), Mobj.nVerts)

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
    patch('Vertices', [x, y] , 'Faces', Mobj.tri, ...
        'Cdata', Mobj.h, 'edgecolor', 'k', 'facecolor', 'interp')
    hold on
    whos Nlist
    plot(x(Nodelist), y(Nodelist), 'ro');
    axis('equal', 'tight')
    title('open boundary nodes');
end

% add to mesh object
npts = numel(Nodelist);
Mobj.nObs = Mobj.nObs + 1;
Mobj.nObcNodes(Mobj.nObs) = npts;
Mobj.obc_nodes(Mobj.nObs,1:npts) = Nodelist;
Mobj.obc_name{Mobj.nObs} = ObcName;
Mobj.obc_type(Mobj.nObs) = ObcType;

if ftbverbose
    fprintf('found %d open boundary nodes', npts)
end

if ftbverbose
    fprintf('\nend   : %s\n', subname)
end

