function write_WRF_forcing(WRF, filename)
% Write data out to WRF format netCDF forcing file.
%
% write_WRF_forcing(WRF, fileprefix)
%
% DESCRIPTION:
%   Takes the given regularly gridded forcing data and writes out to a WRF
%   format netCDF file.
%
% INPUT:
%   WRF - struct with the following fields:
%       lon   : longitude, rectangular array (see MESHGRID).
%       lat   : latitude, rectangular array (see MESHGRID).
%       time  : Modified Julian Day times.
%       pres  : sea level pressure [Pa]
%       nlwrf : net longwave radiation (upward = negative) [W/m^{2}]
%       nswrf : net shortwave radiation (upward = negative) [W/m^{2}]
%       nshf  : net surface heat flux (ocean losing = negative) [W/m^{2}]
%       u10   : eastward wind velocity [m/s]
%       v10   : northward wind velocity [m/s]
%       prate : precipitation (ocean losing = negative) [m/s]
%       evap  : evaporation (ocean losing = negative) [m/s]
%       Net surface heat flux is defined as the sum of shortwave, longwave,
%       sensible and latent heat fluxes (ocean losing heat = negative).
%   filename - Output netCDF file name.
%
% OUTPUT:
%   WRF format heating netCDF file.
%
% EXAMPLE USAGE:
%   wrf_file = '/path/to/output/casename_wnd.nc';
%   write_WRF_forcing(WRF, wrf_file);
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
%   2013-08-29 - First version based on write_FVCOM_heating.m.
%
%==========================================================================

assert(nargin == 2, 'Incorrect number of arguments')

subname = 'write_WRF_heating';

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

ntimes = numel(WRF.time);
[sgYr, sgMon, sgDay, sgHr, sgMin, sgSec] = mjulian2greg(WRF.time(1));
[egYr, egMon, egDay, egHr, egMin, egSec] = mjulian2greg(WRF.time(end));

x = WRF.lon;
y = WRF.lat;
% Make the range of lon -180 - 180.
if max(x) > 180
    x(x > 180) = x(x > 180) - 360;
end
% I've yet to come across a latitude range that isn't -90 - 90, but just in
% case.
if max(y) > 90
    y(y > 90) = y(y > 90) - 180;
end
[nsouth_north, nwest_east] = size(WRF.lon);

%--------------------------------------------------------------------------
% Create the netCDF header for the FVCOM forcing file
%--------------------------------------------------------------------------

nc = netcdf.create(filename, 'clobber');

netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'type', 'FVCOM METEO FORCING FILE')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'title', [filename, ' forcing'])
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'gauge', 'Met Office Unified Model forcing')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'history', ['File created on ', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ' with write_WRF_forcing.m from the MATLAB fvcom-toolbox (https://github.com/pwcazenave/fvcom-toolbox)'])
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'source', 'wrf grid (structured) surface forcing')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'START_DATE', datestr(datenum(sgYr, sgMon, sgDay, sgHr, sgMin, sgSec), 'yyyy-mm-dd HH:MM:SS'))
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'END_DATE', datestr(datenum(egYr, egMon, egDay, egHr, egMin, egSec), 'yyyy-mm-dd HH:MM:SS'))
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'file', filename)
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'institution', 'Plymouth Marine Laboratory')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'Conventions', 'CF-1.0')

% Dimensions
sn_dimid=netcdf.defDim(nc, 'south_north', nsouth_north);
we_dimid=netcdf.defDim(nc, 'west_east', nwest_east);
bu_dimid=netcdf.defDim(nc, 'bottom_up', 1); % not sure if this is necessary...
datestrlen_dimid=netcdf.defDim(nc, 'DateStrLen', 19);
time_dimid=netcdf.defDim(nc, 'time', netcdf.getConstant('NC_UNLIMITED'));

% Space variables
y_varid = netcdf.defVar(nc, 'XLAT', 'NC_FLOAT', [we_dimid, sn_dimid]);
netcdf.putAtt(nc, y_varid, 'long_name', 'latitude');
netcdf.putAtt(nc, y_varid, 'description', 'LATITUDE, SOUTH IS NEGATIVE');
netcdf.putAtt(nc, y_varid, 'units', 'degrees_north');
netcdf.putAtt(nc, y_varid, 'type', 'data');

x_varid=netcdf.defVar(nc, 'XLONG','NC_FLOAT', [we_dimid, sn_dimid]);
netcdf.putAtt(nc, x_varid, 'long_name', 'longitude');
netcdf.putAtt(nc, x_varid, 'description', 'LONGITUDE, WEST IS NEGATIVE');
netcdf.putAtt(nc, x_varid, 'units', 'degrees_east');
netcdf.putAtt(nc, x_varid, 'type', 'data');

% Time variables
times_varid=netcdf.defVar(nc, 'Times', 'NC_CHAR', [datestrlen_dimid, time_dimid]);
netcdf.putAtt(nc, times_varid, 'description', 'GMT time');
% netcdf.putAtt(nc, times_varid, 'format', 'String: Calendar Time');
% netcdf.putAtt(nc, times_varid, 'time_zone', 'UTC');

nswrf_varid = netcdf.defVar(nc, 'Shortwave', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, nswrf_varid, 'long_name', 'Shortwave, upward is negative');
netcdf.putAtt(nc, nswrf_varid, 'units', 'W m-2');
netcdf.putAtt(nc, nswrf_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, nswrf_varid, 'coordinates', 'lat lon');
netcdf.putAtt(nc, nswrf_varid, 'type', 'data');

nshf_varid = netcdf.defVar(nc, 'Shortwave', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, nshf_varid, 'long_name', 'Sum of shortwave, longwave, sensible and latent heat fluxes, ocean lose heat is negative');
netcdf.putAtt(nc, nshf_varid, 'units', 'W m-2');
netcdf.putAtt(nc, nshf_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, nshf_varid, 'coordinates', 'lat lon');
netcdf.putAtt(nc, nshf_varid, 'type', 'data');

u10_varid = netcdf.defVar(nc, 'U10', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, u10_varid, 'long_name', 'Eastward Wind Velocity');
netcdf.putAtt(nc, u10_varid, 'description', 'U at 10 M');
netcdf.putAtt(nc, u10_varid, 'units', 'm s-1');
netcdf.putAtt(nc, u10_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, u10_varid, 'type', 'data');

v10_varid = netcdf.defVar(nc, 'V10', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, v10_varid, 'long_name', 'Northward Wind Velocity');
netcdf.putAtt(nc, v10_varid, 'description', 'V at 10 M');
netcdf.putAtt(nc, v10_varid, 'units', 'm s-1');
netcdf.putAtt(nc, v10_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, v10_varid, 'type', 'data');

prate_varid = netcdf.defVar(nc, 'Precipitation', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, prate_varid, 'long_name', 'Precipitation');
netcdf.putAtt(nc, prate_varid, 'description', 'Precipitation, ocean lose water is negative');
netcdf.putAtt(nc, prate_varid, 'units', 'm s-1');
netcdf.putAtt(nc, prate_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, prate_varid, 'coordinates', 'lat lon');
netcdf.putAtt(nc, prate_varid, 'type', 'data');

evap_varid = netcdf.defVar(nc, 'Evaporation', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, evap_varid, 'long_name', 'Evaporation');
netcdf.putAtt(nc, evap_varid, 'description', 'Evaporation, ocean lose water is negative');
netcdf.putAtt(nc, evap_varid, 'units', 'm s-1');
netcdf.putAtt(nc, evap_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, evap_varid, 'coordinates', 'lat lon');
netcdf.putAtt(nc, evap_varid, 'type', 'data');

pres_varid = netcdf.defVar(nc, 'air_pressure', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, pres_varid, 'long_name', 'Air Pressure');
netcdf.putAtt(nc, pres_varid, 'description', 'Sea surface airpressure');
netcdf.putAtt(nc, pres_varid, 'units', 'Pa');
netcdf.putAtt(nc, pres_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, pres_varid, 'coordinates', 'lat lon');
netcdf.putAtt(nc, pres_varid, 'type', 'data');

% End definitions
netcdf.endDef(nc);

%--------------------------------------------------------------------------
% Put the data in the netCDF file.
%--------------------------------------------------------------------------

% Build the Times string and output to netCDF.
nStringOut = char();
for tt = 1:ntimes
    [nYr, nMon, nDay, nHour, nMin, nSec] = mjulian2greg(data.time(tt));
    nDate = [nYr, nMon, nDay, nHour, nMin, nSec];
    nStringOut = [nStringOut, sprintf('%04i/%02i/%02i %02i:%02i:%02i', nDate)];
end
netcdf.putVar(nc, times_varid, nStringOut);
% And the rest...
netcdf.putVar(nc, x_varid, x);
netcdf.putVar(nc, y_varid, y);
netcdf.putVar(nc, nswrf_varid, WRF.data.nswrf);
netcdf.putVar(nc, nshf_varid, WRF.data.nshf);
netcdf.putVar(nc, u10_varid, WRF.data.u10);
netcdf.putVar(nc, v10_varid, WRF.data.v10);
netcdf.putVar(nc, prate_varid, WRF.data.prate);
netcdf.putVar(nc, evap_varid, WRF.data.evap);
netcdf.putVar(nc, dlwrf_varid, WRF.data.dlwrf.node);
netcdf.putVar(nc, nswrf_varid, WRF.data.dswrf.node);
netcdf.putVar(nc, pres_varid, WRF.data.pres);

% Close the netCDF file(s)
netcdf.close(nc);

if ftbverbose
    fprintf('end   : %s \n', subname)
end
