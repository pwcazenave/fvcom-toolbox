function [Mobj]  = add_river_nodes(Mobj,RiverName)

% Add a set of river nodes comprising a single river to Mesh structure  
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

% use ginput2 (which allows zooming and plots points as they are clicked) to let
% user select the boundary points
fprintf('click the river nodes in the model\n')
fprintf('left click:  zoom\n');
fprintf('delete:  delete last point\n');
fprintf('drag mouse: pan\n');
fprintf('right click: select point\n');
fprintf('enter (return): finished selecting points\n');
[xselect] = ginput2(true,'r+');


[npts,jnk] = size(xselect);

if(npts == 0)
	fprintf('you didn''t select any points')
	fprintf(['end   : ' subname '\n'])
	return
end;
fprintf('you selected %d points\n',npts)

% snap to the closest vertices
for i=1:npts
	[ipt(i),dist] = find_nearest_pt(xselect(i,1),xselect(i,2),Mobj);
end;

% replot domain with snapped vertices
plot(x(ipt),y(ipt),'ro');

% add to mesh object
Mobj.nRivers = Mobj.nRivers + 1;
Mobj.nRivNodes(Mobj.nRivers) = npts;
Mobj.riv_nodes(Mobj.nRivers,1:npts) = ipt;
Mobj.riv_name{Mobj.nRivers} = RiverName;


if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end;

