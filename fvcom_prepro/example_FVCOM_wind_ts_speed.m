function example_FVCOM_wind_ts_speed

% example file for FVCOM, time-varying/spatially constant wind forcing as speed
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
%
% Revision history
%   
%==============================================================================
warning off
subname = 'example_FVCOM_wind_ts';
fprintf('\n')
fprintf(['begin : ' subname '\n'])

%------------------------------------------------------------------------------
% create a dataset
%------------------------------------------------------------------------------

% make a timeframe
% tbeg = greg2julian(2009,1,1,0,0,0)-2400000.5;
% tend = greg2julian(2010,1,1,0,0,0)-2400000.5;
tbeg = 0;
tend = 31;
time = tbeg:1:tend;
nTimes = numel(time);

% make up a fake time varying wind in m/s at 10-m above the water surface
u10 = 10*ones(nTimes,1); 
v10 = zeros(nTimes,1); 

%------------------------------------------------------------------------------
% write output to time-varying, spatially constant FVCOM wind file 
%------------------------------------------------------------------------------
fvcom_forcing_file = 'tst_wind.nc'; 
nc = netcdf(fvcom_forcing_file, 'clobber');            
nc.references = 'http://fvcom.smast.umassd.edu'; 
nc.source = 'single-point time-dependent surface forcing'; 
nc.institution = 'School for Marine Science and Technology' ;
nc.history = 'generated using the fvcom-toolbox';


  
% dimensions
nc('time') = 0;

% time vars
nc{'time'} = ncfloat('time');
nc{'time'}.long_name = 'time';
nc{'time'}.units = 'days since 1858-11-17 00:00:00';
nc{'time'}.format = 'modified julian day (MJD)';
nc{'time'}.time_zone = 'UTC';
  
nc{'Itime'} = ncint('time');
nc{'Itime'}.units = 'days since 1858-11-17 00:00:00';
nc{'Itime'}.format = 'modified julian day (MJD)';
nc{'Itime'}.time_zone = 'UTC';

nc{'Itime2'} = ncint('time');
nc{'Itime2'}.units = 'msec since 00:00:00';
nc{'Itime2'}.time_zone = 'UTC';


nc{'U10'} = ncfloat('time');
nc{'U10'}.long_name = 'Eastward Wind Velocity';
nc{'U10'}.standard_name = 'Wind Velocity';
nc{'U10'}.units = 'm/s';
nc{'U10'}.type = 'data';

nc{'V10'} = ncfloat('time');
nc{'V10'}.long_name = 'Northward Wind Velocity';
nc{'V10'}.standard_name = 'Wind Velocity';
nc{'V10'}.units = 'm/s';
nc{'V10'}.type = 'data';

% dump time
nc{'time'}(1:nTimes) = time; 
nc{'Itime'}(1:nTimes) = floor(time); 
nc{'Itime2'}(1:nTimes) = mod(time,1)*24*3600*1000.;

nc{'U10'}(1:nTimes) = u10;  
nc{'V10'}(1:nTimes) = v10; 

ierr = close(nc);

fprintf(['end   : ' subname '\n'])



