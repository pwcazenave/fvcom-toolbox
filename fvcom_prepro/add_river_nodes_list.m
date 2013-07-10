function [Mobj]  = add_river_nodes(Mobj,Nlist,RiverName)

% Add a set of river nodes comprising a single river to Mesh structure  
% Using a set of user-defined nodes
%
% [Mobj] = add_river_nodes(Mobj)
%
% DESCRIPTION:
%    Select using ginput the set of nodes comprising a river
%
% INPUT
%    Mobj = Matlab mesh object
%    RiverName = Name of the River
%
% OUTPUT:
%    Mobj = Matlab mesh object with an additional river nodelist
%
% EXAMPLE USAGE
%    Mobj = add_river_nodes(Mobj,'Potomac')
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Note:
%    Uses ginput2 which allows zooming before selecting points and displays
%    clicked points realtime
%
% Revision history
%   
%==============================================================================
subname = 'add_river_nodes';
global ftbverbose
if(ftbverbose)
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

%------------------------------------------------------------------------------
% Get a unique list and make sure they are in the range of node numbers 
%------------------------------------------------------------------------------
Nlist = unique(Nlist);

if(max(Nlist) > Mobj.nVerts);
  fprintf('your river node number(s) exceed the total number of nodes in the domain\n');
  fprintf('stop screwing around\n');
  error('stopping...\n')
end;

%------------------------------------------------------------------------------
% Plot the mesh 
%------------------------------------------------------------------------------

if(lower(Mobj.nativeCoords(1:3)) == 'car')
        x = Mobj.x;
        y = Mobj.y;
else
        x = Mobj.lon;
        y = Mobj.lat;
end;

figure
patch('Vertices',[x,y],'Faces',Mobj.tri,...
        'Cdata',Mobj.h,'edgecolor','k','facecolor','interp');
hold on;

plot(x(Nlist),y(Nlist),'ro');
title('river nodes');

% add to mesh object
npts = numel(Nlist);
Mobj.nRivers = Mobj.nRivers + 1;
Mobj.nRivNodes(Mobj.nRivers) = npts;
Mobj.riv_nodes(Mobj.nRivers,1:npts) = Nlist;
Mobj.riv_name{Mobj.nRivers} = RiverName;


if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end;

