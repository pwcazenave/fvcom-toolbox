function [d] = ST_phi2d(phi)
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

%subname = 'ST_phid2';  
%fprintf('\n')
%fprintf(['begin : ' subname '\n'])

%------------------------------------------------------------------------------
% calculate d and convert to m
%------------------------------------------------------------------------------

d = 2^(-phi);
d = d*.001;

%fprintf(['end   : ' subname '\n'])
