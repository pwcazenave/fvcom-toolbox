function example_FVCOM_wind_ts_speed
% Example forcing file for FVCOM, time-varying/spatially constant wind forcing
% as speed.
%
% function example_FVCOM_wind_ts_speed
%
% DESCRIPTION:
%    Write a time-varying, spatially constant wind file
%    This is TEMPLATE program to be used as an example
%    Do not commit your user-defined programs to the repository
%    It requires USER Modification to work 
%
% INPUT
%
% OUTPUT:
%    NetCDF WindFile
%
% EXAMPLE USAGE
%    example_FVCOM_wind_ts_speed
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2016-02-18 Updated to use the relevant functions from the toolbox.
%
%=============================================================================

subname = 'example_FVCOM_wind_ts';
global ftbverbose
ftbverbose = true;
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

%-----------------------------------------------------------------------------
% read a grid in
%-----------------------------------------------------------------------------

Mobj = read_sms_mesh('2dm', './samples/tst.2dm');

%-----------------------------------------------------------------------------
% create a dataset
%-----------------------------------------------------------------------------

% make a timeframe
tbeg = greg2mjulian(2009,1,1,0,0,0);
tend = greg2mjulian(2010,1,1,0,0,0);
data.time = tbeg:1/2:tend;
nTimes = numel(data.time);

% make up a fake time varying wind in m/s at 10-m above the water surface
data.u10.data = 10.*ones(nTimes, Mobj.nElems);
data.v10.data = 5.*ones(nTimes, Mobj.nElems);

%-----------------------------------------------------------------------------
% write output to time-varying, spatially constant FVCOM wind file 
%-----------------------------------------------------------------------------

write_FVCOM_forcing(Mobj, 'tst', data, 'test forcing', '3.2.1')

if ftbverbose
    fprintf('end   : %s\n', subname)
end