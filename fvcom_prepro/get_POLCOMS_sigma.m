function pc = get_POLCOMS_sigma(pc,inputConf)
% Extract bathymetry and sigma layer information from NOC Operation Tide
% Surge model.
%
% function pc = get_POLCOMS_sigma(pc,inputConf)
%
% DESCRIPTION:
%   Extract bathymetry and sigma layer information from NOC Operation Tide
%   Surge model.
%
% INPUT:
%    inputConf   = MATLAB structure which must contain: 
%                   - inputConf.polcoms_ts - location of NOC Operational
%                   Model output containing 4D variables of temperature
%                   (tem) and salinity (sal). They should have dimensions
%                   (x, y, sigma, time).
%                   - inputConf.polcoms_z - location of NOC Operational
%                   Model output containing 4D variables of bathymetry
%                   (XXX) and sigma layer thickness (XXX).
%                   - inputConf.startDate - start date and time for FVCOM
%                   model run
%                   - inputConf.endDate - end date and time for FVCOM
%                   model run
%    pc          = MATLAB structure which must contain:
% 
% OUTPUT:
%    pc     = MATLAB structure to which three new fields (called temperature,
%           salinity and ts_time) have been added.
%
% EXAMPLE USAGE
%    pc = get_POLCOMS_sigma(pc,inputConf)
%
% Author(s):
%    Karen Thurston (National Oceanography Centre, Liverpool)
%
% KJT Revision history:
%    2013-02-05 First version, adapted from Laurent Amoudry's
%    'readparameters.m' and 'plot_bathymetry.m' files.
%
%==========================================================================

subname = 'get_POLCOMS_sigma';

global ftbverbose;
if ftbverbose
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end
%

% Coefficients for the sigma level transformation
hc = 150;  % limit of near-surface high-resolution region
theta = 8;
b = 0.05;
s = 1;  % leave this as 1

% Extract the bathymetry ('pdepth' is cell thickness, 'depth' is cell
% centre depth).
[l,m,n,x,y,D]=readparameters(fullfile(inputConf.polcoms_z,'parameters.dat'),...
    fullfile(inputConf.polcoms_z,'bathymetry.dat'));

% Mask the land points
temp = D==0;
D(temp)=NaN;

% Calculate the sigma levels
[sigma,pc.depth.data,zreg] = scoord3d(D,n,hc,theta,b,s);

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end