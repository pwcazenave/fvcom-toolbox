function [u,v] = naut_2_cart(speed,degT)
% convert wind from nautical (speed,direction) to Cartesian (u,v)
%
% function [u,v] = naut_2_cart(speed,degT)
%
% DESCRIPTION:
%  convert wind from nautical (speed,direction) to Cartesian (u,v)
%
% INPUT:
%   speed:  Wind speed in m/s
%   degT:   Wind direction in DEGREES in Nautical sense (=0, North Wind)
%
% OUTPUT:
%    u:     x-direction (East) wind component in m/s
%    v:     y-direction (North) wind component in m/s
%
% EXAMPLE USAGE
%    [u,v] = naut_2_cart(10,0)
%
% Author(s):  
%   
%
% Revision history
%   
%==============================================================================

deg2rad = pi/180;

angle = (360.0 - degT) + 90.0;
angle = angle*deg2rad;
u = -cos(angle).* speed;
v = -sin(angle).* speed;
