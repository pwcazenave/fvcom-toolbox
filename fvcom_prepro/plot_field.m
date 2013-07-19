function plot_field(Mobj,PlotField,varargin) 

% Plot the mesh, user defined field, open boundary nodes, and river points 
%
% function plot_field(Mobj,field,varargin)
%
% DESCRIPTION:
%    Plot a field from the Mesh structure
%
% INPUT 
%   Mobj                    = matlab mesh object 
%   PlotField               = vertex-based field to plot
%   [optional] coordinate   = coordinate system 
%                             'cartesian'(default)
%                             'spherical'  
%   [optional] showgrid     = show the grid
%                             [true ; false (default)] 
%   [optional] withextra    = display river nodes and obc nodes
%                             [true ; false (default)]
%
% OUTPUT:
%    Figure Plot to Screen
%
% EXAMPLE USAGE
%    plot_field(Mobj,Mobj.h,'coordinate','spherical')
%    plot_field(Mobj,Mobj.ts,'showgrid',true)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

subname = 'plot_field';
global ftbverbose
if(ftbverbose);
  fprintf('\n'); fprintf(['begin : ' subname '\n'])
end;

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
PlotCartesian = 'cartesian';
ShowGrid = false;
HaveTitle = false;
PlotExtra = false;

nArgs = length(varargin);
if(mod(nArgs,2) ~= 0)
	error('incorrect usage of plot_mesh')
end;


for i=1:2:length(varargin)-1
	keyword  = lower(varargin{i});
	if( ~ischar(keyword) )
		error('incorrect usage of plot_mesh')
	end;
	
	switch(keyword(1:3))
	
	case'coo'
		coord = lower(varargin{i+1});
		if(coord(1:3) == 'car')
			PlotCartesian = true;
		else
			PlotCartesian = false;
		end;
	case'sho'
		showg = lower(varargin{i+1});
		if(showg)
			ShowGrid = true;
		else
			ShowGrid = false;
		end;
	case'tit'
		MyTitle = varargin{i+1};
		HaveTitle = true;
	case'wit'
		showg = lower(varargin{i+1});
		if(showg)
			PlotExtra = true;
		else
			PlotExtra = false;
		end;
	otherwise
		error(['Can''t understand value for:' keyword]);
	end; %switch keyword
end;

%------------------------------------------------------------------------------
% Plot the mesh and bathymetry
%------------------------------------------------------------------------------

if(ShowGrid)
	edgecolor = 'k';
else
	edgecolor = 'interp';
end;

field = PlotField;

if(PlotCartesian)
	if(~Mobj.have_xy)
		error('no (x,y) coordinates available for Mesh structure')
	end;
	x = Mobj.x;
	y = Mobj.y;
else
	if(~Mobj.have_lonlat)
		error('no (lon,lat) coordinates available for Mesh structure')
	end;
	x = Mobj.lon;
	y = Mobj.lat;
end;

figure
patch('Vertices',[x,y],'Faces',Mobj.tri,...
      	'Cdata',field,'edgecolor',edgecolor,'facecolor','interp');
colorbar
hold on;

if(HaveTitle)
	title(MyTitle)
end;

%------------------------------------------------------------------------------
% Plot the mesh and bathymetry
%------------------------------------------------------------------------------
if(PlotExtra)
	for i=1:Mobj.nRivers
		nodes = Mobj.riv_nodes(i,1:Mobj.nRivNodes(i));
		plot(x(nodes),y(nodes),'go','MarkerFaceColor','g');
	end;
	for i=1:Mobj.nObs
		nodes = Mobj.obc_nodes(i,1:Mobj.nObcNodes(i));
		plot(x(nodes),y(nodes),'ks','MarkerFaceColor','k');
	end;
	for i=1:Mobj.nSponge
		nodes = Mobj.sponge_nodes(i,1:Mobj.nSpongeNodes(i));
		plot(x(nodes),y(nodes),'rd','MarkerSize',5,'MarkerFaceColor','k');
	end;
end;


if(ftbverbose);
  fprintf(['end   : ' subname '\n'])
end;


