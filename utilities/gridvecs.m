function gridvecs(gridfile,x1,x2,y1,y2,ds,thresh,fname)

%clear all; close all;
%gridfile = '~gcowles/Data/NOAA_MITSG_TKE/mtm1b.nc';
%x1 = 836100;
%x2 = 935370;
%y1 = -186230;
%y2 = -125810;
%ds = 2000;
%fname = 'mtm1b_2km_nansound_vecs.dat';

% Find elements closest to a grid of points for cleaning up vector plots 
%
% function [] = gridvecs(fname,x2,y1,y2,ds)
%
% DESCRIPTION:
% Find elements closest to a grid of points for cleaning up vector plots 
%
% INPUT:
%    gridfile:  FVCOM netcdf output file containing x,y   
%    fname:  list of node numbers 
%    x1:  left side of bounding box
%    x2:  right side of bounding box
%    y1:  bottom of bounding box
%    y2:  top of bounding box
%    ds:  grid lengthscale
%    thresh: threshold distance in units of x.  If nearest FVCOM point to the 
%            vector grid is greater than this distance, no point is assigned.
%    
% OUTPUT:
%    fname: file containing indices of node numbers where vectors should be plot
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% References
%==============================================================================


% open grid file and read data 
if(~exist(gridfile))
  error('grid file does not exist')
end;

nc = netcdf(gridfile);
x  = nc{'x'}(:);
y  = nc{'y'}(:);
nVerts = numel(x);
fprintf('number of nodes in the mesh %d\n',nVerts);
nc = close(nc);

% create the grid
[X,Y] = meshgrid(x1:ds:x2,y1:ds:y2);
[il,jl] = size(X);
node = zeros(il*jl,1);
cnt   = 0;
for i=1:il
for j=1:jl
   xloc = X(i,j); 
   yloc = Y(i,j);
   dist = sqrt(   (x-xloc).^2 + (y-yloc).^2);  
   [mind,imin] = min(dist);
   if(mind < thresh);  
     cnt = cnt + 1;
     node(cnt) = imin;
   end;
end;
end;
node = node(1:cnt);  

% dump the list
fid = fopen(fname,'w');
for i=1:cnt
  fprintf(fid,'%d\n',node(i));
end;
fclose(fid);
