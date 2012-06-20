function [phi] = ST_d2phi(d)
% Convert sediment diameter (m) to phi 
%
% function [phi] = ST_d2phi(d)
%
% DESCRIPTION:
%    Convert a sediment grain size in m to phi 
%
% INPUT:
%    d: sediment grain size in m
%
% OUTPUT:
%    phi: sediment grain size in phi scale
%
% EXAMPLE USAGE
%    phi = ST_d2phi(.001)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

%subname = 'ST_d2phi';  
%fprintf('\n')
%fprintf(['begin : ' subname '\n'])

%------------------------------------------------------------------------------
% calculate phi 
%------------------------------------------------------------------------------
phi = -log2(d*1000);

%fprintf(['end   : ' subname '\n'])
