function example_FVCOM_tsobc(basename, time, nSiglay)
% example file for dumping a file to force temperature and salinity at the
% open boundary.
%
% function example_FVCOM_tsobc(basename, time, nSiglay)
%
% DESCRIPTION:
%    Setup a sample FVCOM hydrographic open boundary forcing file
%
% INPUT
%    Model case name
%    Time (Modified Julian Days)
%    Number of sigma layers
%
% OUTPUT:
%    FVCOM hydrographic open boundary file
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2012-06-15 Added support for native MATLAB NetCDF routines. Requires
%    MATLAB 2010a or higher.
%    2012-07-16 Removed hard-coded nSiglay and nSiglev and instead moved to
%    arguments list.
%    2016-02-18 Updated the output to netCDF to use the relevant function.
%
%==============================================================================

subname = 'example_FVCOM_tsobc';
global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

fvcom_bathy = [basename, '_dep.dat'];
fvcom_obc   = [basename, '_obc.dat'];
tsOBCFile = [basename, '_tsobc.nc'];

nTimes = numel(time);

% Set siglev/siglay
nSiglev = nSiglay + 1;
inc = 1./real(nSiglay);
siglev = 0:-inc:-1;
for i = 1:nSiglay
    siglay(i) = mean(siglev(i:i+1));
end

% Initialize temperature/salinity arrays
% temp = zeros(nObc,nSiglay,nTimes);
% salt = zeros(nObc,nSiglay,nTimes);
% set variable temperature and salinity
% for i=1:nTimes
%     obc_temp(i) = 18. + 2.*real(i-1)/nTimes;
%     obc_salt(i) = 30. - 5.*real(i-1)/nTimes;
% end

% set a constant temperature and salinity
obc_temp = 13;
obc_salt = 35;

%--------------------------------------------------------------
% dump to netcdf file
%--------------------------------------------------------------

Mobj.siglev = siglev;
Mobj.siglay = siglay;
write_FVCOM_tsobc(basename,time,nSiglay,obc_temp,obc_salt,Mobj)

if ftbverbose
    fprintf('end   : %s\n', subname);
end
