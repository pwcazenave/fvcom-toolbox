function [Mobj]  = add_sponge_nodes(Mobj,SpongeName,SpongeRadius,SpongeCoeff)

% Add a set of sponge nodes comprising a single sponge layer to Mesh structure  
%
% [Mobj] = add_sponge_nodes(Mobj)
%
% DESCRIPTION:
%    Select using ginput the set of nodes comprising a sponge layer
%
% INPUT
%    Mobj = Matlab mesh object
%    SpongeName = Name of the Sponge Layer
%    SpongeRadius = Radius of influence of the Sponge Layer 
%    SpongeCoeff  = Sponge damping coefficient
%
% OUTPUT:
%    Mobj = Matlab mesh object with an additional sponge nodelist
%
% EXAMPLE USAGE
%    Mobj = add_sponge_nodes(Mobj,'Sponge1',10000,.0001)
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
subname = 'add_sponge_nodes';
fprintf('\n')
fprintf(['begin : ' subname '\n'])


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
axis('equal','tight')
hold on;

% use ginput2 (which allows zooming and plots points as they are clicked) to let
% user select the boundary points
[xselect] = ginput2(true,'k+')


npts = size(xselect,1);

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
Mobj.nSponge = Mobj.nSponge + 1;
Mobj.nSpongeNodes(Mobj.nSponge) = npts;
Mobj.sponge_nodes(Mobj.nSponge,1:npts) = ipt;
Mobj.sponge_name{Mobj.nSponge} = SpongeName;
Mobj.sponge_rad(Mobj.nSponge) = SpongeRadius;
Mobj.sponge_fac(Mobj.nSponge) = SpongeCoeff;


fprintf(['end   : ' subname '\n'])

