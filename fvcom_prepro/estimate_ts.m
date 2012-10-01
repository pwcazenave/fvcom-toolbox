function [Mobj] = estimate_ts(Mobj,u,zeta)

% Estimate time step at each node
%
% [Mobj] = function estimate_ts(Mobj)
%
% DESCRIPTION:
%    Calculate barotropic time step
%
% INPUT
%    Mobj = matlab mesh object
%
% OUTPUT:
%    Mobj = matlab structure containing mesh time step
%
% EXAMPLE USAGE
%    Mobj = estimate_ts(Mobj)
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%    2012-07-14 Add great circle approximation if only provided with
%    latitude and longitudes. Also add arguments to the function to define
%    current velocity and tidal amplitudes.
%
%==============================================================================

subname = 'estimate_ts';
global ftbverbose
if(ftbverbose)
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end;

%------------------------------------------------------------------------------
% Set constants
%------------------------------------------------------------------------------

g    = 9.81; %gravitational acceleration
%u    = 3.0;  %u-velocity
%zeta = 11.0;  %tide amp

if(~Mobj.have_bath)
    error('can''t estimate the time step without bathymetry')
end;

%------------------------------------------------------------------------------
% Compute the time step estimate
%------------------------------------------------------------------------------
if Mobj.have_xy
    x = Mobj.x;
    y = Mobj.y;
else
    % Will convert to metres when calculating element edge length
    x = Mobj.lon;
    y = Mobj.lat;
end
h = max(Mobj.h,0);
tri = Mobj.tri;
nVerts = Mobj.nVerts;
nElems = Mobj.nElems;

ts = ones(nVerts,1)*1e9;
lside = zeros(nVerts,1);
for i=1:nElems
    n1 = tri(i,1);
    n2 = tri(i,2);
    n3 = tri(i,3);
    nds = [n1 n2 n3];
    % Check whether we have x and y values and use great circle
    % approximations if we don't.
    if Mobj.have_xy
        lside(i) = sqrt( (x(n1)-x(n2))^2 + (y(n1)-y(n2))^2);
    else
        lside(i) = haversine(x(n1),y(n1),x(n2),y(n2));
    end
    dpth  = max(h(nds))+zeta;
    dpth  = max(dpth,1);
    ts(nds) = min(ts(nds),lside(i)/(sqrt(g*dpth) + u));
end;
if(ftbverbose); fprintf('minimum time step: %f seconds\n',min(ts)); end;
Mobj.ts = ts;
Mobj.have_ts = true;

if(ftbverbose)
    fprintf(['end   : ' subname '\n'])
end

function [km]=haversine(lat1,lon1,lat2,lon2)
% Haversine function to calculate first order distance measurement. Assumes
% spherical Earth surface. Lifted from:
%
% http://www.mathworks.com/matlabcentral/fileexchange/27785
%
R = 6371000;                    % Earth's mean radius in metres
delta_lat = lat2 - lat1;        % difference in latitude
delta_lon = lon2 - lon1;        % difference in longitude
a = sin(delta_lat/2)^2 + cos(lat1) * cos(lat2) * ...
    sin(delta_lon/2)^2;
c = 2 * atan2(sqrt(a), sqrt(1-a));
km = R * c;                     % distance in metres

