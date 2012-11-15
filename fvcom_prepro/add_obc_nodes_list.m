function [Mobj]  = add_obc_nodes_list(Mobj,Nlist,ObcName,ObcType) 

% Add a set of obc nodes comprising a single obc boundary to Mesh structure  
% Using a list of nodes
%
% [Mobj] = add_obc_nodes_list(Mobj,Nlist,ObcName,ObcType)
%
% DESCRIPTION:
%    Select using ginput the set of nodes comprising an obc
%
% INPUT
%    Mobj = Matlab mesh object
%    Nlist = List of nodes
%    ObcName = Name of the Open Boundary
%    ObcType = FVCOM Flag for OBC Type
%
% OUTPUT:
%    Mobj = Matlab mesh object with an additional obc nodelist
%
% EXAMPLE USAGE
%    Mobj = add_obc_nodes_list(Mobj,Nlist,'OpenOcean')
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
%
% Revision history:
%   
%==========================================================================
subname = 'add_obc_nodes';
global ftbverbose
if(ftbverbose)
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end


%--------------------------------------------------------------------------
% Get a unique list and make sure they are in the range of node numbers 
%--------------------------------------------------------------------------
Nlist = unique(Nlist);

if(max(Nlist) > Mobj.nVerts);
  fprintf('your open boundary node number exceed the total number of nodes in the domain\n');
  error('stopping...\n')
end

%--------------------------------------------------------------------------
% Plot the mesh 
%--------------------------------------------------------------------------

if strcmpi(Mobj.nativeCoords(1:3), 'car')
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
whos Nlist
plot(x(Nlist),y(Nlist),'ro');
axis('equal','tight')
title('open boundary nodes');

% add to mesh object
npts = numel(Nlist);
Mobj.nObs = Mobj.nObs + 1;
Mobj.nObcNodes(Mobj.nObs) = npts; 
Mobj.obc_nodes(Mobj.nObs,1:npts) = Nlist;
Mobj.obc_name{Mobj.nObs} = ObcName;
Mobj.obc_type(Mobj.nObs) = ObcType;


if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end

