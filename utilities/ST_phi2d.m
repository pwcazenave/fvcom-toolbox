function d = ST_phi2d(phi)
% Convert sediment grain size from phi to d (m)
%
% function [d] = ST_phi2d(phi)
%
% DESCRIPTION:
%    Convert a sediment grain size from phi to d (m)
%
% INPUT:
%    phi: sediment grain size in phi scale
%
% OUTPUT:
%    d: sediment grain size in m
%
% EXAMPLE USAGE
%    d = ST_phi2d(2.5)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

global ftbverbose
[~, subname] = fileparts(mfilename('fullpath'));
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

%------------------------------------------------------------------------------
% calculate d and convert to m
%------------------------------------------------------------------------------

d = 2^(-phi);
d = d*.001;

if ftbverbose
    fprintf('end   : %s\n', subname)
end