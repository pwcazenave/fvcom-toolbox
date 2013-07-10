function [cell]  = inCell(Mobj,xpt,ypt)

% Find cell in mesh Mobj containing point (xpt,ypt)
%
% function [cell]  = inCell(Mobj,xpt,ypt)
%
% DESCRIPTION:
%   Find cell in mesh Mobj containing point (xpt,ypt)
%
% INPUT
%    Mobj = Matlab mesh object
%    xpt = x position of point
%    ypt = y position of point
%
% OUTPUT:
%    cell = cell number (=0 if not found in mesh) 
%
% EXAMPLE USAGE
%    cell = inCell(Mobj,5.,6.)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
tri    = Mobj.tri;
nElems = Mobj.nElems;
x      = Mobj.x;
y      = Mobj.y;
cell   = 0;
xc     = Mobj.xc;
yc     = Mobj.yc;
dist   = sqrt(  (xc-xpt).^2 + (yc-ypt).^2);

cell = 0;
%try nearest
[rmin,imin] = min(dist); 
if(isintriangle(x(tri(imin,1:3)),y(tri(imin,1:3)),xpt,ypt));  
  cell = imin; 
  return; 
end;

%sort from min distance to max and search along sort
[distsort,ind] = sort(dist,1,'ascend');
for i=1:nElems
   if(isintriangle(x(tri(ind(i),1:3)),y(tri(ind(i),1:3)),xpt,ypt));  
     cell = ind(i); 
     return
   end;
end;
  
% didn't find cell, just return nearest cell
%cell = ind(1);
