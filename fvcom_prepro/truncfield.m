function [field]  = truncfield(fieldin,Mobj,scale) 

% Smooth a vertex-based field using averages  
%
% [field] = function smoothfield(fieldin,Mobj,SmoothFactor,nLoops,SmoothPts)  
%
% DESCRIPTION:
%    Smooth a vertex based field 
%
% INPUT
%    Mobj         = Matlab mesh object
%    fieldin      = vertex based field
%    scale        = scale in grid coordinates over which to calc trunc value 
%
% OUTPUT:
%    field = filtered, vertex-based field
%
% EXAMPLE USAGE
%    Mobj.h = smoothfield(Mobj.h,Mobj,1000)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

%------------------------------------------------------------------------------
% Parse input
%------------------------------------------------------------------------------

if(exist('fieldin')*exist('Mobj')*exist('scale') == 0)
	error('arguments to truncfield are missing')
end;



%------------------------------------------------------------------------------
% Smoothing Loops
%------------------------------------------------------------------------------
m = Mobj.nVerts;
for i=1:m
  radlist(1:m,1) = sqrt((Mobj.x(1:m)-Mobj.x(i)).^2 + (Mobj.y(1:m)-Mobj.y(i)).^2);
  ipts = find(radlist < scale);
  fieldmin = min(abs(fieldin(ipts)));
  field(i) = fieldmin;
end;

