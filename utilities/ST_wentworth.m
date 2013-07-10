function [Sclass] = wentworth(phi)
% Report wentworth class of a particular grain size phi 
%
% function wentworth(phi) 
%
% DESCRIPTION:
%    Report the wentworth of a grain size phi
%
% INPUT:
%    phi: sediment grain size in phi scale
%
% OUTPUT:
%    prints wentworth class to screen
%
% EXAMPLE USAGE
%    wentworth(2.5) 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

ClassNames = {'boulder','cobble','pebble','granule','very coarse sand', ...
   'coarse sand','medium sand','fine sand','very fine sand','coarse silt', ...
   'medium silt','fine silt','very fine silt','coarse clay','medium clay','fine clay'}; 
ClassLbound = [-1e6,-8,-6,-2,-1,0,1,2,3,4,5,6,7,8,9,10]; 
pts = find(phi-ClassLbound > 0);
ClassIndex = pts(end); 
Sclass = char(ClassNames{ClassIndex});
%fprintf('class of phi = %f is: %s\n',phi,class);
