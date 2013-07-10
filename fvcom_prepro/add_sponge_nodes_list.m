function [Mobj]  = add_sponge_nodes_list(Mobj,SpongeList,SpongeName,SpongeRadius,SpongeCoeff,plotFig)

% Add a set of sponge nodes comprising a single sponge layer to Mesh structure  
%
% [Mobj] = add_sponge_nodes(Mobj)
%
% DESCRIPTION:
%    Select using ginput the set of nodes comprising a sponge layer
%
% INPUT
%    Mobj = Matlab mesh object
%    SpongeList = List of nodes to which to create a Sponge Layer
%    SpongeName = Name of the Sponge Layer
%    SpongeRadius = Radius of influence of the Sponge Layer 
%    SpongeCoeff  = Sponge damping coefficient
%    plotFig = [optional] show a figure of the mesh (1 = yes)
%
% OUTPUT:
%    Mobj = Matlab mesh object with an additional sponge nodelist
%
% EXAMPLE USAGE
%    Mobj = add_sponge_nodes(Mobj,'Sponge1',10000,.0001)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Karen Thurston (National Oceanography Centre, Liverpool)
%
% Revision history
%    Modifed from add_sponge_nodes to read in nodes from a supplied list.
%    2012-11-26 Add ability to turn off the figures.
%    2013-01-18 Added support for variable sponge radius
%   
%==============================================================================
subname = 'add_sponge_nodes';

global ftbverbose
if(ftbverbose)
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end

% Do we want a figure showing how we're getting along?
if nargin == 5
    plotFig = 0;
end

%------------------------------------------------------------------------------
% Plot the mesh 
%------------------------------------------------------------------------------

if plotFig == 1
    if(lower(Mobj.nativeCoords(1:3)) == 'car')
        x = Mobj.x;
        y = Mobj.y;
    else
        x = Mobj.lon;
        y = Mobj.lat;
    end

    figure
    patch('Vertices',[x,y],'Faces',Mobj.tri,...
            'Cdata',Mobj.h,'edgecolor','k','facecolor','interp');
    hold on;
    plot(x(SpongeList),y(SpongeList),'wx')
    axis('equal','tight')
end

npts = length(SpongeList);

if(npts == 0)
	fprintf('No points in given list')
	fprintf(['end   : ' subname '\n'])
	return
end
if(ftbverbose)
    fprintf('%d points provided\n',npts)
end

% add to mesh object
Mobj.nSponge = Mobj.nSponge + 1;
Mobj.nSpongeNodes(Mobj.nSponge) = npts;
Mobj.sponge_nodes(Mobj.nSponge,1:npts) = SpongeList;
Mobj.sponge_name{Mobj.nSponge} = SpongeName;
Mobj.sponge_fac(Mobj.nSponge) = SpongeCoeff;

if max(size(SpongeRadius))==1   % if you have a constant sponge radius
    Mobj.sponge_rad(Mobj.nSponge) = SpongeRadius;
else    % if you have a variable sponge radius
    Mobj.sponge_rad(Mobj.nSponge,1:npts) = SpongeRadius;
end

if(ftbverbose)
    fprintf(['end   : ' subname '\n'])
end

