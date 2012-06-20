function [Mobj]  = add_sponge_nodes_list(Mobj,SpongeList,SpongeName,SpongeRadius,SpongeCoeff)

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
%
% Revision history
%    Modifed from add_sponge_nodes to read in nodes from a supplied list. 
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
hold on;
plot(x(SpongeList),y(SpongeList),'wx')
axis('equal','tight')

npts = size(SpongeList,2);

if(npts == 0)
	fprintf('No points in given list')
	fprintf(['end   : ' subname '\n'])
	return
end;
fprintf('%d points provided\n',npts)

% add to mesh object
Mobj.nSponge = Mobj.nSponge + 1;
Mobj.nSpongeNodes(Mobj.nSponge) = npts;
Mobj.sponge_nodes(Mobj.nSponge,1:npts) = SpongeList;
Mobj.sponge_name{Mobj.nSponge} = SpongeName;
Mobj.sponge_rad(Mobj.nSponge) = SpongeRadius;
Mobj.sponge_fac(Mobj.nSponge) = SpongeCoeff;


fprintf(['end   : ' subname '\n'])

