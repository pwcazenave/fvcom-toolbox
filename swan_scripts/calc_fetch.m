function [dir,fetch] = calc_fetch(nodefile,botfile,xp,yp,zadd);
% estimate wave fetch 
%
% function [dir,fetch] = calc_fetch(nodefile,botfile,xp,yp,zadd);
%
% DESCRIPTION:
%   given a SWAN mesh and bathymetry and observation point, estimate fetch for  
%   range of wind directions
%
% INPUT 
%   nodefile = SWAN unstructured node file
%   botfile  = SWAN unstructure bathymetry file (positive down)
%   xp       = observation point x-coordinate
%   yp       = observation point y-coordinate
%   zadd     = [optional] add free surface height to account for tides (m)
%
% OUTPUT:
%    dir = direction of wind in Cartesian coordinates
%    fetch = fetch for that direction in (m)
%
%    Note:  if dir=45, this will provide the fetch for a SW wind
%
% EXAMPLE USAGE
%   [dir,fetch] = calc_fetch('skg4.3.node','skg4.3.bot',5.3870e5,5.3506e6,-2);
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

%clear all; close all;
%nodefile = 'skg4.3.node';
%elefile = 'skg4.3.ele';
%botfile = 'skg4.3.bot';
%xp = 5.3870e+05;
%yp = 5.3506e+06;
%zeta_max = 2.5;
%zeta_min = 2.5;

% set number of theta-bins
nbins = 32; 

delta_zeta = 0.0;
if(exist('zadd'))
  delta_zeta = zadd;
end;

% warning - program is terribly unrobust
% to improve robustness, search along rays until the point is found
% inside a boundary or dry triangle

% read the node file
[cnt,x,y,nmark] = textread(nodefile,'%d %f %f %d\n','headerlines',1);

% read the depth file
[h] = textread(botfile,'%f\n','headerlines',0);

% add/subtract free surface amplitude
h = h + delta_zeta;

% mark all nodes that are solid wall, open boundary, or dry
drynodes = find(h < 0);
nmark(drynodes) = 3;
solidpts = find(nmark > 0);
xdry = x(solidpts);
ydry = y(solidpts);

% distance vector from all solidpts to point of interest
dist = sqrt(   (xdry-xp).^2 + (ydry-yp).^2);
ang  = atan2( (ydry-yp), (xdry-xp) );


% find the wave direction corresponding to the angle between the boundary point and POI
plot(dist.*cos(ang),dist.*sin(ang),'r+'); hold on;

% sort the points in order of increasing angle 
ang = ang + pi; %convert from obs-shore to shore-obs angle

% for each bin find the minimum distance
dtheta = 2*pi/nbins;
theta_bound = 0-dtheta/2:dtheta:2*pi+dtheta/2;
theta       = 0:dtheta:2*pi;
mind        = zeros(numel(theta),1);
for i=1:numel(theta)
  pts = find( (theta_bound(i) < ang) & (ang <= theta_bound(i+1)));
  if(numel(pts) > 0)
    mind(i) = min(dist(pts)); 
  end;
end;

% make sure there are no zeros 
for i=2:numel(theta)-1
  if(mind(i) == 0.); mind(i) = max(mind(i-1),mind(i+1)); end;
end;
if(mind(1) == 0.); mind(1) = mind(2); end;
if(mind(end) == 0.); mind(end) = mind(end-1); end;

% smooth to account for diffraction
for i=2:numel(theta)-1
  mind(i) = .5*(mind(i) + max(mind(i-1),mind(i+1))); 
end;
mind(1) = .5*(mind(1) + mind(2));
mind(end) = .5*(mind(end) + mind(end-1));

% force periodicity
mind(end) = mind(1);
  
plot(0,0,'g+','MarkerSize',5); hold on;
plot(mind.*cos(theta+pi)',mind.*sin(theta+pi)');

% plot a polar   
figure
plot(theta*180/pi,mind);
xlabel('wind direction (Cartesian sense)');
ylabel('fetch (m)');
 
% set return values
dir = theta*180/pi;
fetch = mind;


