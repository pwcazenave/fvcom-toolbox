function example_FVCOM_wind_ts

% example file for FVCOM, time-varying/spatially constant wind forcing as stress
%
% function example_FVCOM_wind_ts
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
%    example_FVCOM_wind_ts
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
warning off
subname = 'example_FVCOM_wind_ts';
global ftbverbose;
if(ftbverbose);
fprintf('\n')
fprintf(['begin : ' subname '\n'])
end;

%------------------------------------------------------------------------------
% create a dataset
%------------------------------------------------------------------------------

% make a timeframe
% tbeg = greg2julian(2009,1,1,0,0,0)-2400000.5;
% tend = greg2julian(2010,1,1,0,0,0)-2400000.5;
tbeg = 0;
tend = 31;
time = tbeg:(1./24.):tend;
nTimes = prod(size(time));

% make up a fake time varying wind
taux = 0.25*(cos( ((time-time(1))*2*pi)/7)) + .15*(cos( ((time-time(1))*2*pi)/360));
tauy = 0.25*(sin( ((time-time(1))*2*pi)/7)) + .15*(sin( ((time-time(1))*2*pi)/360));

% plot the wind
subplot(2,1,1)
plot(time-time(1),taux,'k'); 
subplot(2,1,2)
plot(time-time(1),tauy,'r');

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


nc{'uwind_stress'} = ncfloat('time');
nc{'uwind_stress'}.long_name = 'Eastward Wind Stress'; 
nc{'uwind_stress'}.standard_name = 'Wind Stress'; 
nc{'uwind_stress'}.units = 'Pa';
nc{'uwind_stress'}.type = 'data';

nc{'vwind_stress'} = ncfloat('time');
nc{'vwind_stress'}.long_name = 'Northward Wind Stress'; 
nc{'vwind_stress'}.standard_name = 'Wind Stress'; 
nc{'vwind_stress'}.units = 'Pa';
nc{'vwind_stress'}.type = 'data';

% dump time
nc{'time'}(1:nTimes) = time; 
nc{'Itime'}(1:nTimes) = floor(time); 
nc{'Itime2'}(1:nTimes) = mod(time,1)*24*3600*1000.;

nc{'uwind_stress'}(1:nTimes) = taux; 
nc{'vwind_stress'}(1:nTimes) = tauy; 

ierr = close(nc);

if(ftbverbose);
fprintf(['end   : ' subname '\n'])
end;



