function [Mobj]  = add_obc_nodes_graphic(Mobj,ObcName,ObcType)

% Add a set of obc nodes comprising a single obc boundary to Mesh structure  
% By clicking on points on the screen
%
% [Mobj] = add_obc_nodes(Mobj)
%
% DESCRIPTION:
%    Select using ginput the set of nodes comprising an obc
%
% INPUT
%    Mobj = Matlab mesh object
%    ObcName = Name of the Open Boundary
%    ObcType = FVCOM Flag for OBC Type
%
% OUTPUT:
%    Mobj = Matlab mesh object with an additional obc nodelist
%
% EXAMPLE USAGE
%    Mobj = add_obc_nodes(Mobj,'OpenOcean')
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Note:
%    Uses ginput2 which allows zoom/pan before selecting points and displays
%    clicked points realtime
%
% Revision history
%   
%==============================================================================
subname = 'add_obc_nodes';
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
[xselect] = ginput2(true,'k+')


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
Mobj.nObs = Mobj.nObs + 1;
Mobj.nObcNodes(Mobj.nObs) = npts;
Mobj.obc_nodes(Mobj.nObs,1:npts) = ipt;
Mobj.obc_name{Mobj.nObs} = ObcName;
Mobj.obc_type(Mobj.nObs) = ObcType;


if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end;

