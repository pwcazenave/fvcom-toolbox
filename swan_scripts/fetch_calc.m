function [fetch_struct] = fetch_calc(fvcom_mesh,fvcom_bath,lon,lat,zeta_min,zeta_max,nZeta,nTheta,dry_thresh,darc) 
% close all; clear all;
% fvcom_mesh = 'skg4.3_grd.dat';
% fvcom_bath = 'skg4.3_dep.dat';
% lon = -122.4721646 ;
% lat = 48.3372476; 
% zeta_min = -3;
% zeta_max = 3;
% nZeta = 30;
% nTheta = 16;
% dry_thresh = 0.1;
% darc = 125;
% Calculate fetch as a function of free surface height (positive up) and orientation 
%
% function [fetch_struct] = fetch_calc(fvcom_mesh,lon,lat,zeta_min,zeta_max,dzeta,dtheta,dry_thresh) 
%
% DESCRIPTION:
%   Calculate fetch as a function of free surface height (positive up) and orientation 
%
% INPUT 
%   fvcom_mesh  = FVCOM 3.x grid file 
%   fvcom_bath  = FVCOM 3.x bathymetry file 
%   lon         = longitude of point
%   lat         = latitude of point
%   zeta_min    = minimum value of free surface height (m)
%   zeta_max    = maximum value of free surface height (m)
%   nZeta       = number of zeta partitions 
%   nTheta      = number of theta partitions 
%   dry_thresh  = threshold for a point to be dry
%   darc        = arc spacing in meters along which fetch is searched
%
% OUTPUT:
%   fetch_struct = fetch structure containing fetch as a function of angle and free surface 
%                  height.  Angle is defined using Cartesian Wind Stress.  Thus the fetch 
%                  corresponding to an angle of pi/2 (90 degrees) is the
%                  fetch corresponding to a South wind (wind from South
%                  with Northward wind stress)
%
% EXAMPLE USAGE
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

%-------------------------------------------------
% load the mesh metrics
%-------------------------------------------------

% read the Mesh from an FVCOM grid file
Mobj = read_fvcom_mesh(fvcom_mesh); 
x = Mobj.x;
y = Mobj.y;

% read the bathymetry from an FVCOM grid file
Mobj.h = read_fvcom_bath(fvcom_bath); 
h = Mobj.h;

% set cell center coordinates
Mobj.xc = nodes2elems(x,Mobj);
Mobj.yc = nodes2elems(y,Mobj);

%-------------------------------------------------
% set grids in zeta and theta 
%-------------------------------------------------
theta = -pi:2*pi/nTheta:pi; nTheta=numel(theta);
zeta  = zeta_min:(zeta_max-zeta_min)/nZeta:zeta_max;
if(zeta_min==zeta_max)
  zeta = zeta_min;
  nZeta = 1;
else
  nZeta = nZeta+1;
end;

%-------------------------------------------------
% project observation point (lon,lat) => (x,y) 
%-------------------------------------------------
[xobs,yobs] = my_project(lon,lat,'forward'); 


%-------------------------------------------------
% compute elevation (positive up) of observation point   
%-------------------------------------------------
zobs = -griddata(x,y,h,xobs,yobs); 

%-------------------------------------------------
% loop over angles and depths, compute fetch 
%-------------------------------------------------

radline = 0:darc:max(max(x)-min(x),max(y)-min(y)); nRad = numel(radline);
fetch = zeros(nTheta,nZeta);
hc = nodes2elems(Mobj.h,Mobj);
for i=1:nTheta
  fprintf('searching theta %f\n',theta(i)*180/pi);
  xline = xobs+radline*cos(theta(i)+pi);
  yline = yobs+radline*sin(theta(i)+pi);
  k = 1;
  for j=1:nZeta
    elev = zeta(j); depth = hc+elev;
    found = 0;
    k = 1;
    while ~found
      cell = inCell(Mobj,xline(k),yline(k));
      if(cell == 0)| (depth(cell) < dry_thresh);
        fetch(i,j) = radline(k);
        %fprintf('found at k = %d\n',k)
        found = true;
      end
      k = k + 1;
    end;
  end;
end;
           

%-------------------------------------------------
% save to a structure 
%-------------------------------------------------

fetch_struct.xobs  = xobs;
fetch_struct.yobs  = yobs;
fetch_struct.zobs  = zobs;
fetch_struct.theta = theta;
fetch_struct.zeta  = zeta;
fetch_struct.fetch = fetch;
fetch_struct.xmesh = x;
fetch_struct.ymesh = y;
fetch_struct.hmesh = h;


