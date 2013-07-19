function [Mobj] = gen_unstructured_mesh(lx,ly,res,cornermod) 

% function gen_unstructured_mesh(lx,ly,res,cornermod)
%
% DESCRIPTION:
%    Generate an unstructured mesh in a rectangular domain of size dist_x*dist_y 
%    having a resolution of res.  
%
% INPUT
%    lx        :  domain length in x-direction
%    ly        :  domain length in y-direction
%    res       :  target triangle edge length in same units as lx,ly 
%    cornermod :  ensure no cell shares a boundary on two of the rectangles edges 
%
%   
% OUTPUT:
%    Preprocessing toolbox mesh object
%
% EXAMPLE USAGE
%   gen_unstructured_mesh(1000,100,10,false); 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% TODO
%   Implement the edgeswap to satisfy cornermod=true
%
%
%  Dependency - Matlab mesh2d - at least v24
% http://www.mathworks.com/matlabcentral/fileexchange/25555-mesh2d-automatic-mesh-generation
%    
%==============================================================================

subname = 'gen_unstructured_mesh';
fprintf('\n')
fprintf(['begin : ' subname '\n'])

userproject = false;

% test vals
% lx = 1000;
% ly = 500;
% cornermod = true;
% res = 50;
% close all; warning off;

% generate mesh
node = [0 0; lx 0; lx ly ; 0 ly];
edge = [1 2 ; 2 3 ; 3 4; 4 1]; 
hdata.hmax = res; 
options.output = true;
[p,t] = mesh2d(node,edge,hdata,options);

% generate connectivity
[e,te,et2,bnd] = connectivity(p,t);


% option to eliminate corner elements TODO
% if(cornermod)
%     
%   for pt=1:1
%     xpt = node(pt,1); ypt = node(pt,2);
%     dist = sqrt( (xpt-p(:,1)).^2 + (ypt-p(:,2)).^2);
%     [imin,ipt] = min(dist);
%   end; 
% 
% 
% end;

%------------------------------------------------------------------------------
% Create a blank mesh object
%------------------------------------------------------------------------------
Mobj = make_blank_mesh();
coordinate = 'cartesian';




% ensure all elements are CCW
area = triarea(p,t);
fprintf('min area %f max area %f\n',min(area),max(area));

%------------------------------------------------------------------------------
% Transfer to Mesh structure
%------------------------------------------------------------------------------

% set dimensions
nElems = numel(t(:,1));
nVerts = numel(p(:,1));
fprintf('# of elements in the mesh %d\n',nElems);
fprintf('# of verts    in the mesh %d\n',nVerts);

% allocate memory to hold mesh and auxiliary structures
tri = zeros(nElems,3);
x   = zeros(nVerts,1);
y   = zeros(nVerts,1);
h   = zeros(nVerts,1);
lon = zeros(nVerts,1);
lat = zeros(nVerts,1);
ts  = zeros(nVerts,1);


Mobj.nVerts  = nVerts;
Mobj.nElems  = nElems;
Mobj.nativeCoords = coordinate;

Mobj.have_lonlat  = false;
Mobj.have_xy      = true;
Mobj.have_bath    = false;

Mobj.x            = p(:,1);
Mobj.y            = p(:,2);
Mobj.ts           = ts;
Mobj.lon          = lon;
Mobj.lat          = lat;
Mobj.h            = h;
Mobj.tri          = t;


fprintf(['end   : ' subname '\n'])
%write_dtank('junk.dtascii',Mobj);


