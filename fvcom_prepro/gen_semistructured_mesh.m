function [Mobj] = gen_semistructured_mesh(lx,ly,nx,ny,cornermod) 

% Generate a semistructured mesh 
%
% [Mobj] = gen_semistructured_mesh(varargin)
%
% DESCRIPTION:
%    Read SMS 2dm file and bathymetry file 
%    Store in a matlab mesh object 
%
% INPUT 
%   lx  = length of the domain in meters in the x-direction
%   ly  = length of the domain in meters in the y-direction
%   nx  = number of edges along the x-direction 
%   ny  = number of edges in the y-direction
%   cornermod = modify edges at the corners so that one cell does not have an 
%   edge on two sides (true,false)
%
% OUTPUT:
%    Mobj = matlab structure containing mesh data
%
% EXAMPLE USAGE
%    Mobj = gen_semistructured_mesh(10000.,1000.,500,50)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

subname = 'gen_semistructured_mesh';
fprintf('\n')
fprintf(['begin : ' subname '\n'])

userproject = false;

%------------------------------------------------------------------------------
% Create a blank mesh object
%------------------------------------------------------------------------------
Mobj = make_blank_mesh();
coordinate = 'cartesian';

%------------------------------------------------------------------------------
% set grid resolution
%------------------------------------------------------------------------------
dx = lx/real(nx);
dy = ly/real(ny);
fprintf('making a semistructured mesh of dimensions: %f X %f\n',lx,ly)
fprintf('resolution in x %f and y %f\n',dx,dy);


%------------------------------------------------------------------------------
% build the grid
%------------------------------------------------------------------------------

il = nx+1;
jl = ny+1;

% set dimensions
nElems = 2*nx*ny;
nVerts = il*jl;

% allocate memory to hold mesh and auxiliary structures
tri = zeros(nElems,3);
x   = zeros(nVerts,1);
y   = zeros(nVerts,1);
h   = zeros(nVerts,1);
lon = zeros(nVerts,1);
lat = zeros(nVerts,1);
ts  = zeros(nVerts,1);


% build the mesh points
nn = 0;
for i=1:il
  for j=1:jl
    nn = nn + 1;
    x(nn) = (i-1)*dx;
    y(nn) = (j-1)*dy;
  end;
end;

% set the connectivity
nn = 1;
for i=1:nx
	for j=1:ny
		n2 = nn + 1;
		tri(nn,1) = j + (i-1)*jl;
		tri(nn,2) = j + i*jl;
		tri(nn,3) = j+1 + (i-1)*jl;
		tri(n2,1) = j + i*jl;
		tri(n2,2) = j+1 + i*jl;
		tri(n2,3) = j+1 + (i-1)*jl;
		nn = nn + 2;
	end;
end;

		
% modify the connectivity on the corners to deal with the constraint that an element cannot
% have both an open and a solid boundary
if(cornermod)
	tri(1,1) = 1;
	tri(1,2) = jl+2;
	tri(1,3) = 2;
	tri(2,1) = 1;
	tri(2,2) = jl+1;
	tri(2,3) = jl+2;

	tri(2*nx*ny-1,1) = il*jl-jl-1;
	tri(2*nx*ny-1,2) = il*jl-1;
	tri(2*nx*ny-1,3) = il*jl;
	tri(2*nx*ny  ,1) = il*jl-jl-1;
	tri(2*nx*ny  ,2) = il*jl;
	tri(2*nx*ny  ,3) = il*jl-jl;
end;


%------------------------------------------------------------------------------
% Transfer to Mesh structure
%------------------------------------------------------------------------------

Mobj.nVerts  = nVerts;
Mobj.nElems  = nElems;
Mobj.nativeCoords = coordinate;

Mobj.have_lonlat  = false;
Mobj.have_xy      = true;
Mobj.have_bath    = false;

Mobj.x            = x;
Mobj.y            = y;
Mobj.ts           = ts;
Mobj.lon          = lon;
Mobj.lat          = lat;
Mobj.h            = h;
Mobj.tri          = tri;


fprintf(['end   : ' subname '\n'])


