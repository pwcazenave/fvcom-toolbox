function [instant_rad, rad] = calc_rad(cloud, lat, jul, secs)
% Calculating instantaneous irradiance from cloud cover.
%
% [instant_rad, rad] = calc_rad(cloud, lat, jul, secs)
%
% DESCRIPTION:
%   Calculate the surface irradiance (short wave) as a function of
%   latitude and date.
%
% INPUT:
%   cloud   - cloud cover in percentage (0-1)
%   lat     - latitude of observations
%   jul     - MATLAB date (see HELP DATENUM)
%   secs    - Seconds in the day
%
% OUTPUT:
%   instant_rad     -
%   rad             -
% 
% EXAMPLE USAGE:
%   [cloud, lat, jul, secs] = deal(0.8, 42, floor(datenum(T_TIME(12))), 43200)
%
% NOTES:
%   This routine is taken from GOTM-ERSEM.
%
% Author(s):
%   Ricard Torres (Plymouth Marine Laboratory)
%
% Revision history:
%   2014-01-14 Added to the FVCOM toolbox.
%   What units!!!?

subname = 'calc_rad';

global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

cycle = 365.0;
pi = 3.141592653589793;
rfactor = pi/180.0;
rconstant = 1368.0;
dq = 149.6;
aa = 149.499;
eps = 0.01675;
sk = 1368.0;
noon = 12.0;

% Parameter values
a = [0.4, 0.517, 0.474, 0.421, 0.38, 0.35, 0.304, 0.23, 0.106, 0.134];
b = [0.386, 0.317, 0.381, 0.413, 0.468, 0.457, 0.438, 0.384, 0.285, 0.295];
nmth = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 30, 30];
lmth = [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 30, 30];
cnt = -1;
rcnt = -1;

% Obtain actual day
[yy, mm, dd, hh, min, ss] = datevec(jul);

% Leap year?
leap_switch = 2000 - yy;
leap_switch = mod(leap_switch, 4);

switch leap_switch
    case 0
        day = 0;
        for j = 1:mm
            day = day + lmth(j);
        end
        day = day + dd;
    otherwise
        day = 0;
        for j = 1:mm
            day = day + nmth(j);
        end
        day = day + dd;
end
atime = day - 1 + secs / (3600.0 * 24.0);

latirad = lat * rfactor;

corr = 0.09 * cloud * 8.0;

declination = -0.406 * cos(2.0* pi * floor(atime) / cycle);
daylength = acos(-tan(declination) * tan(latirad));

rad = -1.0 * ((sin(declination) * sin(latirad) * daylength / 2.0 + ...
    cos(declination) * cos(latirad) * sin(daylength / 2.0)) * ...
    sk / pi * (1.0 - corr));
daylength = daylength / pi * 24.0;

t = (secs - noon * 3600.0) * pi / (12.0 * 3600.0);
t0 = (noon - 0.5 * daylength) * 3600.0;
t1 = (noon + 0.5 * daylength) * 3600.0;
instant_rad = rad * (cos(t / (2 * daylength / 24.0)) + ...
    0.5 * (1 + cos(t / (daylength / 24.0)))) * 2 * pi / (4 + pi);
if (secs < t0) || (secs > t1)
    instant_rad = 0.0;
end

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end

return
