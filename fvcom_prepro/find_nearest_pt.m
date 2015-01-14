function [Point,Distance] = find_nearest_pt(xloc,yloc,Mobj) 
% Find nearest point in Mesh structure to (x,y) 
%
% function [Point,Distance] = find_nearest_pt(xloc,yloc,Mobj) 
%
% DESCRIPTION:
%    Find nearest point to (xloc,yloc) in the domain of Mobj
%    using native coordinates of Mobj
%
% INPUT:
%   xloc   = x location of point (in native Mobj coordinates)
%   yloc   = y location of point (in native Mobj coordinates)
%   Mobj   = Mesh object with the following fields:
%               - nativeCoords = grid type (cartesian or spherical)
%               - x, y and/or lon, lat = coordinates (dependent on
%               nativeCoords).
%
% OUTPUT:
%   Point = index of nearest vertex in the mesh
%   Distance = Distance from x,y to Point in Mobj native coordinates
%
% EXAMPLE USAGE
%    [Point,Distance] = find_nearest_point(50.1,100.2,Mobj)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2015-01-14 Tidy up the code a bit and add extra information to the
%    help.
%
%==============================================================================

global ftbverbose
subname = 'find_nearest_pt';
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
if ~exist('xloc', 'var') || ~exist('yloc', 'var') || ~exist('Mobj', 'var')
	error('arguments to %s are missing', subname)
end

%------------------------------------------------------------------------------
% Set native coordinates
%------------------------------------------------------------------------------
if strcmpi(Mobj.nativeCoords, 'cartesian')
	x = Mobj.x;
	y = Mobj.y;
elseif strcmpi(Mobj.nativeCoords, 'spherical')
	x = Mobj.lon;
	y = Mobj.lat;
else
    error('Unrecognised coordinate type.')
end

%------------------------------------------------------------------------------
% Find the nearest point
%------------------------------------------------------------------------------
[Distance, Point] = min(sqrt((xloc - x).^2 + (yloc - y).^2));

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end


