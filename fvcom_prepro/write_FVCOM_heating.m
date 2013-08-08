function write_FVCOM_heating(Mobj, fileprefix, data)
% Write data out to FVCOM NetCDF heating file (for HEATING_CALCULATED).
%
% write_FVCOM_heating(Mobj, fileprefix, data)
%
% DESCRIPTION:
%   Takes the given interpolated data (e.g. from grid2fvcom) and writes out
%   to a NetCDF file.
%
% INPUT:
%   Mobj - MATLAB mesh object containing fields:
%       tri - triangulation table for the unstructured grid
%       nVerts - number of grid vertices (nodes)
%       nElems - number of grid elements
%       nativeCoords - model coordinate type ('cartesian' or 'spherical')
%       x, y or lon, lat - node positions (depending on nativeCoords value)
%   fileprefix - Output NetCDF file prefix (plus path) will be
%       fileprefix_{wnd,hfx,evap}.nc if fver is '3.1.0', otherwise output
%       files will be fileprefix_wnd.nc.
%   data - Struct of the data to be written out.
%
% The fields in data may be called any of:
%     - 'slp'               - sea level pressure
%     - 'rhum'              - relative humidity
%     - 'dlwrf'             - downward longwave radiation
%     - 'dswrf'             - downward shortwave radiation
%     - 'air'               - air temperature
%     - 'lon'               - longitude (vector)
%     - 'lat'               - latitude (vector)
%     - 'x'                 - eastings (vector)
%     - 'y'                 - northings (vector)
%
% OUTPUT:
%   FVCOM heating NetCDF file.
%
% EXAMPLE USAGE:
%   heatBase = '/path/to/output/casename_wnd.nc';
%   write_FVCOM_heating(Mobj, heatBase, data);
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
%   2013-07-18 - First version based on write_FVCOM_forcing.m.
%
%==========================================================================

assert(nargin == 3, 'Incorrect number of arguments')

subname = 'write_FVCOM_heating';

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

tri = Mobj.tri;
nNodes = Mobj.nVerts;
nElems = Mobj.nElems;
ntimes = numel(data.time);

if strcmpi(Mobj.nativeCoords, 'cartesian')
    x = Mobj.x;
    y = Mobj.y;
else
    x = Mobj.lon;
    y = Mobj.lat;
    % Make the range of lon 0-360
    x(x < 0) = x(x < 0) + 360;
end
% Create a string for each variable's coordinate attribute
coordString = sprintf('FVCOM %s coordinates', Mobj.nativeCoords);

% Create element vertices positions
xc = nodes2elems(x, Mobj);
yc = nodes2elems(y, Mobj);

%--------------------------------------------------------------------------
% Create the NetCDF header for the FVCOM forcing file
%--------------------------------------------------------------------------

nc = netcdf.create([fileprefix, '_hfx.nc'], 'clobber');

netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'title', 'FVCOM Forcing File')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'institution', 'Plymouth Marine Laboratory')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'source', 'FVCOM grid (unstructured) surface forcing')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'history', ['File created on ', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ' with write_FVCOM_forcing.m from the MATLAB fvcom-toolbox (https://github.com/pwcazenave/fvcom-toolbox)'])
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'references', 'http://fvcom.smast.umassd.edu, http://codfish.smast.umassd.edu')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'Conventions', 'CF-1.0')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'CoordinateSystem', Mobj.nativeCoords)
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'CoordinateProjection', 'init=epsg:4326') % WGS84?

% Dimensions
nele_dimid=netcdf.defDim(nc, 'nele', nElems);
node_dimid=netcdf.defDim(nc, 'node', nNodes);
three_dimid=netcdf.defDim(nc, 'three', 3);
time_dimid=netcdf.defDim(nc, 'time', netcdf.getConstant('NC_UNLIMITED'));
datestrlen_dimid=netcdf.defDim(nc, 'DateStrLen', 26);

% Space variables
if strcmpi(Mobj.nativeCoords, 'cartesian')
    x_varid=netcdf.defVar(nc, 'x', 'NC_FLOAT', node_dimid);
    netcdf.putAtt(nc, x_varid, 'long_name', 'nodal x-coordinate');
    netcdf.putAtt(nc, x_varid, 'units', 'meters');

    y_varid=netcdf.defVar(nc, 'y', 'NC_FLOAT', node_dimid);
    netcdf.putAtt(nc, y_varid, 'long_name', 'nodal y-coordinate');
    netcdf.putAtt(nc, y_varid, 'units', 'meters');

    xc_varid=netcdf.defVar(nc, 'xc','NC_FLOAT', nele_dimid);
    netcdf.putAtt(nc, xc_varid, 'long_name', 'zonal x-coordinate');
    netcdf.putAtt(nc, xc_varid, 'units', 'meters');

    yc_varid=netcdf.defVar(nc, 'yc','NC_FLOAT', nele_dimid);
    netcdf.putAtt(nc, yc_varid, 'long_name', 'zonal y-coordinate');
    netcdf.putAtt(nc, yc_varid, 'units', 'meters');

elseif strcmpi(Mobj.nativeCoords, 'spherical')
    x_varid=netcdf.defVar(nc, 'lon','NC_FLOAT', node_dimid);
    netcdf.putAtt(nc, x_varid, 'long_name', 'nodal longitude');
    netcdf.putAtt(nc, x_varid, 'units', 'degrees_east');

    y_varid = netcdf.defVar(nc, 'lat', 'NC_FLOAT', node_dimid);
    netcdf.putAtt(nc, y_varid, 'long_name', 'nodal latitude');
    netcdf.putAtt(nc, y_varid, 'units', 'degrees_north');

    xc_varid = netcdf.defVar(nc, 'lonc', 'NC_FLOAT', nele_dimid);
    netcdf.putAtt(nc, xc_varid, 'long_name', 'zonal longitude');
    netcdf.putAtt(nc, xc_varid, 'units', 'degrees_east');

    yc_varid = netcdf.defVar(nc, 'latc', 'NC_FLOAT', nele_dimid);
    netcdf.putAtt(nc, yc_varid, 'long_name', 'zonal latitude');
    netcdf.putAtt(nc, yc_varid, 'units', 'degrees_north');   
else
    error('Unknown coordinate type (%s)', Mobj.nativeCoords)
end

nv_varid=netcdf.defVar(nc, 'nv', 'NC_INT', [nele_dimid, three_dimid]);
netcdf.putAtt(nc, nv_varid, 'long_name', 'nodes surrounding element');

% Time variables
time_varid=netcdf.defVar(nc, 'time', 'NC_FLOAT', time_dimid);
netcdf.putAtt(nc, time_varid, 'long_name', 'time');
netcdf.putAtt(nc, time_varid, 'units', 'days since 1858-11-17 00:00:00');
netcdf.putAtt(nc, time_varid, 'format', 'modified julian day (MJD)');
netcdf.putAtt(nc, time_varid, 'time_zone', 'UTC');

itime_varid=netcdf.defVar(nc, 'Itime', 'NC_INT', time_dimid);
netcdf.putAtt(nc, itime_varid, 'units', 'days since 1858-11-17 00:00:00');
netcdf.putAtt(nc, itime_varid, 'format', 'modified julian day (MJD)');
netcdf.putAtt(nc, itime_varid, 'time_zone', 'UTC');

itime2_varid=netcdf.defVar(nc, 'Itime2', 'NC_INT', time_dimid);
netcdf.putAtt(nc, itime2_varid, 'units', 'msec since 00:00:00');
netcdf.putAtt(nc, itime2_varid, 'time_zone', 'UTC');
netcdf.putAtt(nc, itime2_varid, 'long_name', 'time');

times_varid=netcdf.defVar(nc, 'Times', 'NC_CHAR', [datestrlen_dimid, time_dimid]);
netcdf.putAtt(nc, times_varid, 'long_name', 'Calendar Date');
netcdf.putAtt(nc, times_varid, 'format', 'String: Calendar Time');
netcdf.putAtt(nc, times_varid, 'time_zone', 'UTC');

airt_varid = netcdf.defVar(nc, 'air_temperature', 'NC_FLOAT', [node_dimid, time_dimid]);
netcdf.putAtt(nc, airt_varid, 'long_name', 'Surface air temperature');
netcdf.putAtt(nc, airt_varid, 'units', 'Celsius Degree');
netcdf.putAtt(nc, airt_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, airt_varid, 'coordinates', coordString);
netcdf.putAtt(nc, airt_varid, 'type', 'data');

rhum_varid = netcdf.defVar(nc, 'relative_humidity', 'NC_FLOAT', [node_dimid, time_dimid]);
netcdf.putAtt(nc, rhum_varid, 'long_name', 'surface air relative humidity');
netcdf.putAtt(nc, rhum_varid, 'units', 'percentage');
netcdf.putAtt(nc, rhum_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, rhum_varid, 'coordinates', coordString);
netcdf.putAtt(nc, rhum_varid, 'type', 'data');

dlwrf_varid = netcdf.defVar(nc, 'long_wave', 'NC_FLOAT', [node_dimid, time_dimid]);
netcdf.putAtt(nc, dlwrf_varid, 'long_name', 'Downward solar longwave radiation flux');
netcdf.putAtt(nc, dlwrf_varid, 'units', 'Watts meter-2');
netcdf.putAtt(nc, dlwrf_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, dlwrf_varid, 'coordinates', coordString);
netcdf.putAtt(nc, dlwrf_varid, 'type', 'data');

dswrf_varid = netcdf.defVar(nc, 'short_wave', 'NC_FLOAT', [node_dimid, time_dimid]);
netcdf.putAtt(nc, dswrf_varid, 'long_name', 'Downward solar shortwave radiation flux');
netcdf.putAtt(nc, dswrf_varid, 'units', 'Watts meter-2');
netcdf.putAtt(nc, dswrf_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, dswrf_varid, 'coordinates', coordString);
netcdf.putAtt(nc, dswrf_varid, 'type', 'data');

slp_varid = netcdf.defVar(nc, 'air_pressure', 'NC_FLOAT', [node_dimid, time_dimid]);
netcdf.putAtt(nc, slp_varid, 'long_name', 'Surface air pressure');
netcdf.putAtt(nc, slp_varid, 'units', 'Pa');
netcdf.putAtt(nc, slp_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, slp_varid, 'coordinates', coordString);
netcdf.putAtt(nc, slp_varid, 'type', 'data');

% End definitions
netcdf.endDef(nc);

%--------------------------------------------------------------------------
% Put the data in the NetCDF file.
%--------------------------------------------------------------------------
netcdf.putVar(nc, nv_varid, tri');
netcdf.putVar(nc, time_varid, 0, ntimes, data.time);
netcdf.putVar(nc, itime_varid, 0, ntimes, floor(data.time));
netcdf.putVar(nc, itime2_varid, 0, ntimes, mod(data.time, 1) * 24 * 3600 * 1000);
netcdf.putVar(nc, x_varid, x);
netcdf.putVar(nc, y_varid, y);
netcdf.putVar(nc, xc_varid, xc);
netcdf.putVar(nc, yc_varid, yc);
netcdf.putVar(nc, airt_varid, data.air.node);
netcdf.putVar(nc, rhum_varid, data.rhum.node);
netcdf.putVar(nc, dlwrf_varid, data.dlwrf.node);
netcdf.putVar(nc, dswrf_varid, data.dswrf.node);
try % work with both slp and pres data.
    netcdf.putVar(nc, slp_varid, data.slp.node);
catch
    netcdf.putVar(nc, slp_varid, data.pres.node);
end

% Build the Times string and output to NetCDF.
nStringOut = char();
for tt=1:ntimes
    [nYr, nMon, nDay, nHour, nMin, nSec] = mjulian2greg(data.time(tt));
    nDate = [nYr, nMon, nDay, nHour, nMin, nSec];
    nStringOut = [nStringOut, sprintf('%04i/%02i/%02i %02i:%02i:%02i       ', nDate)];
end
netcdf.putVar(nc, times_varid, nStringOut);

% Close the NetCDF file(s)
netcdf.close(nc);

if ftbverbose
    fprintf('end   : %s \n', subname)
end
