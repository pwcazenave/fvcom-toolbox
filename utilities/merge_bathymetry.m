% Blend two sources of bathymetry 
function [ new_bat ] = merge_bathymetry( dist_OB,coarse_bat,fine_bat,d0,d1 )
%
% [ new_bat ] = merge_bathymetry( dist_OB,coarse_bat,fine_bat,d0,d1 )
%
% DESCRIPTION:
%   This script uses a blending function to combine two bathymetry sources
%    over a common spatial extent (High-resolution, unstructured meshes for hydrodynamic models of the Great Barrier Reef, Australia,
%    Estuarine Coastal and Shelf Science 68:36-46 · June 2006
%    DOI: 10.1016/j.ecss.2005.08.017 )

% INPUT [keyword pairs]:
%  dist_OB:     a length metric used to combine the two sources
%  coarse_bat:  Bathymetry data to combine at lower distance 
%  fine_bat:    Bathymetry data to combine at longer distances
%  d0:          is the plateau distance
%  d1:          is the end distance limiting factor
%
%
%  the blending function is a piecewise polynomial function as 
%  f1d= 3.*( (dist_OB-d0)./(d1-d0) ).^2 - 2.*( (dist_OB-d0)./(d1-d0) ).^3;
% % restric closes points
% f1d(dist_OB<=d0)=0;
% f1d(dist_OB>=d1)=1;
% f1d=abs(f1d-1);
% % % Example 
% d0=50;d1=2000;dist_OB = 0:100:12000;
% f1d= 3.*( (dist_OB-d0)./(d1-d0) ).^2 - 2.*( (dist_OB-d0)./(d1-d0) ).^3;
% % restric closes points
% f1d(dist_OB<=d0)=0;
% f1d(dist_OB>=d1)=1;
% plot(dist_OB,abs(f1d-1));
% EXAMPLE USAGE
%    Mobj_all.sub.z_merged = merge_bathymetry( Distance_from_nodestring,Mobj_all.z,Mobj_all.z1,100,2000 );
%
% Author(s):
%    Ricardo Torres (Plymouth Marine Laboratory)
%
% Revision history
%
%   
%==============================================================================



f1d= 3.*( (dist_OB-d0)./(d1-d0) ).^2 - 2.*( (dist_OB-d0)./(d1-d0) ).^3;
% restric closes points
f1d(dist_OB<=d0)=0;
f1d(dist_OB>=d1)=1;
f1d=abs(f1d-1);
new_bat=f1d.*coarse_bat+((1-f1d).*fine_bat);
return

