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
%   Mobj   = Mesh object
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
%
% Revision history
%   
%==============================================================================

subname = 'find_nearest_pt';
%fprintf('\n')
%fprintf(['begin : ' subname '\n'])

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
if(exist('xloc')*exist('yloc')*exist('Mobj') == 0)
	error('arguments to find_nearest_pt are missing')
end;

%------------------------------------------------------------------------------
% Set native coordinates
%------------------------------------------------------------------------------
if(lower(Mobj.nativeCoords(1:3)) == 'car')
	x = Mobj.x;
	y = Mobj.y;
else
	x = Mobj.lon;
	y = Mobj.lat;
end;

%------------------------------------------------------------------------------
% Find the nearest point
%------------------------------------------------------------------------------
radvec = sqrt( (xloc-x).^2 + (yloc-y).^2);
[Distance,Point] = min(radvec);


%fprintf(['end   : ' subname '\n'])


