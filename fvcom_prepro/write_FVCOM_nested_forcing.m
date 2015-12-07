function write_FVCOM_nested_forcing(nest, ncfile, nesttype)
% Creates an FVCOM nesting file.
%
% function write_FVCOM_nested_forcing(nest, ncfile, nesttype)
%
% DESCRIPTION:
%   Uses timeseries data from structured grid already interpolated into
%   FVCOM nodes and elements and generates a netCDF file to drive FVCOM at
%   boundaries
%
% Optionally specify nesting type:
%   1/2: DIRECT/INDIRECT nesting:
%       - Full variables/no surface elevation respectively.
%   3:   RELAXATION nesting:
%       - Nesting with a relaxation method.
%
% INPUT:
%   nest     = struct whose field names are the variable names to be
%              included in netCDF file. Additional required fields are:
%               - time (in Modified Julian Days)
%               - nv (triangulation table)
%               - lon, lat, x, y - node coordinates (spherical and
%               cartesian).
%               - lonc, latc, xc, yc - element coordinates (spherical and
%               cartesian).
%   ncfile   = full path to the nesting file to be created.
%   nesttype = [optional] nesting type (defaults to 1 = direct nesting).
%
% OUTPUT:
%   FVCOM nesting file.
%
% EXAMPLE USAGE:
%   nest.time = time [modified julian days (floats)]
%   nest.temp = Temperature [degrees celsius]
%   nest.salinity = Salinity [PSU]
%   nest.ua = Vertically averaged x velocity [ms^{-1}]
%   nest.va = Vertically averaged y velocity [ms^{-1}]
%   nest.u = Eastward Water Velocity [ms^{-1}]
%   nest.v = Northward Water Velocity [ms^{-1}]
%   nest.hyw = hydro static vertical velocity [ms^{-1}?]
%   nest.weight_cell = weights see manual for explanation [0-1]
%   nest.weight_node = weights see manual for explanation [0-1]
%
%   write_FVCOM_nested_forcing(nest, '/tmp/fvcom_nested.nc', 1)
%
% Author(s):
%   Ricardo Torres (Plymouth Marine Laboratory)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%   Darren Price (CH2MHill)
%   Hakeem Johnson (CH2MHill)
%
% Revision history:
%   2013-06-04 First version based on Riqui's write_nesting_struct_fvcom.m
%   script.
%   2015-02-19 Updated to use either weighted or non-weighted nesting. Also
%   general tidy up.
%   2015-02-24 Add extra time variables, which although not strictly
%   necessary for the model to run, makes sanity checking the time series a
%   lot more straightforward.
%
%==========================================================================

% We need the following variables:
%
% zeta:         Sea surface elevation           [node, time]
% ua:           Vertically averaged x velocity  [node, time]
% va:           Vertically averaged y velocity  [nele, time]
% u:            Eastward Water Velocity         [nele, siglay, time]
% v:            Northward Water Velocity        [nele, siglay, time]
% temp:         Temperature                     [node, siglay, time]
% salinity:     Salinity                        [node, siglay, time]
% hyw:          Hydro static vertical velocity  [node, siglev, time]
% weight_cell:  Weighting for elements          [nele]
% weight_node:  Weighting for nodes             [node]
% Itime:        Days since 1858-11-17 00:00:00  [time]
% Itime2:       msec since 00:00:00             [time]
%
% We include these optional ones for humans:
% time:         Modified Julian Day             [time]
% Times:        Gregorian dates                 [time, datestrlen]

subname = 'write_FVCOM_nested_forcing';

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

if nargin == 2
    nesttype = 1;
elseif nargin < 2 || nargin > 3
    error(['Incorrect input arguments. Supply netCDF file path, ', ...
        'nesting struct and optionally the nesting type (1, 2 or 3).'])
end

% Check we have all the data we need.
required = {'time', 'x', 'y', 'lon', 'lat', 'xc', 'yc', 'lonc', 'latc', ...
    'nv', 'h', 'u', 'v', 'ua', 'va', 'temp', 'salinity', 'hyw', ...
    'weight_cell', 'weight_node'};
fields = fieldnames(nest);
for f = required
    if any(strcmpi(f{1}, {'weight_node', 'weight_cell'})) && nesttype == 3
        assert(any(strcmpi(f, fields)), 'Missing %s struct field', f{1});
    elseif any(strcmpi(f{1}, {'weight_node', 'weight_cell'})) && nesttype ~= 3
        continue
    end
end

[elems, nsiglay, ntimes] = size(nest.u);
nsiglev = nsiglay + 1;
[nodes, ~] = size(nest.zeta);

nc = netcdf.create(ncfile, 'clobber');

% define global attributes
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'type', ...
    'FVCOM nestING TIME SERIES FILE')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'title', ...
    sprintf('FVCOM nestING TYPE %d TIME SERIES data for open boundary', ...
    nesttype))
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'history', ...
    sprintf('File created using %s from the MATLAB fvcom-toolbox', subname))
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'filename', ncfile)
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'Conventions', 'CF-1.0')

% define dimensions
elem_dimid = netcdf.defDim(nc, 'nele', elems);
node_dimid = netcdf.defDim(nc, 'node', nodes);
three_dimid = netcdf.defDim(nc, 'three', 3);
time_dimid = netcdf.defDim(nc, 'time', netcdf.getConstant('NC_UNLIMITED'));
siglay_dimid = netcdf.defDim(nc, 'siglay', nsiglay);
siglev_dimid = netcdf.defDim(nc, 'siglev', nsiglev);
datestrlen_dimid = netcdf.defDim(nc, 'DateStrLen', 26);

% define variables
time_varid = netcdf.defVar(nc, 'time', 'NC_FLOAT', time_dimid);
netcdf.putAtt(nc, time_varid, 'long_name', 'time');
netcdf.putAtt(nc, time_varid, 'units', 'days since 1858-11-17 00:00:00');
netcdf.putAtt(nc, time_varid, 'format', 'modified julian day (MJD)');
netcdf.putAtt(nc, time_varid, 'time_zone', 'UTC');

itime_varid = netcdf.defVar(nc, 'Itime', 'NC_INT', ...
    time_dimid);
netcdf.putAtt(nc, itime_varid, 'units', 'days since 1858-11-17 00:00:00');
netcdf.putAtt(nc, itime_varid, 'format', 'modified julian day (MJD)');
netcdf.putAtt(nc, itime_varid, 'time_zone', 'UTC');

itime2_varid = netcdf.defVar(nc, 'Itime2', 'NC_INT', ...
    time_dimid);
netcdf.putAtt(nc, itime2_varid, 'units', 'msec since 00:00:00');
netcdf.putAtt(nc, itime2_varid, 'time_zone', 'UTC');

Times_varid=netcdf.defVar(nc, 'Times' ,'NC_CHAR', ...
    [datestrlen_dimid, time_dimid]);
netcdf.putAtt(nc, Times_varid, 'time_zone', 'UTC');

x_varid = netcdf.defVar(nc, 'x', 'NC_FLOAT', ...
    node_dimid);
netcdf.putAtt(nc, x_varid, 'units', 'meters');
netcdf.putAtt(nc, x_varid, 'long_name', 'nodal x-coordinate');

y_varid = netcdf.defVar(nc, 'y', 'NC_FLOAT', ...
    node_dimid);
netcdf.putAtt(nc, y_varid, 'units', 'meters');
netcdf.putAtt(nc, y_varid, 'long_name', 'nodal y-coordinate');

xc_varid = netcdf.defVar(nc, 'xc', 'NC_FLOAT', ...
    elem_dimid);
netcdf.putAtt(nc, xc_varid, 'units', 'meters');
netcdf.putAtt(nc, xc_varid, 'long_name', 'zonal x-coordinate');

nv_varid = netcdf.defVar(nc, 'nv', 'NC_INT', ...
    [elem_dimid, three_dimid]);
netcdf.putAtt(nc, xc_varid, 'units', 'no units');
netcdf.putAtt(nc, xc_varid, 'long_name', 'elements nodes indices');

yc_varid = netcdf.defVar(nc, 'yc', 'NC_FLOAT', ...
    elem_dimid);
netcdf.putAtt(nc, yc_varid, 'units', 'meters');
netcdf.putAtt(nc, yc_varid, 'long_name', 'zonal y-coordinate');

lon_varid = netcdf.defVar(nc, 'lon', 'NC_FLOAT', ...
    node_dimid);
netcdf.putAtt(nc, lon_varid, 'units', 'degrees_east');
netcdf.putAtt(nc, lon_varid, 'standard_name', 'longitude');
netcdf.putAtt(nc, lon_varid, 'long_name', 'nodal longitude');

lat_varid = netcdf.defVar(nc, 'lat', 'NC_FLOAT', ...
    node_dimid);
netcdf.putAtt(nc, lat_varid, 'units', 'degrees_north');
netcdf.putAtt(nc, lat_varid, 'standard_name', 'latitude');
netcdf.putAtt(nc, lat_varid, 'long_name', 'nodal latitude');

lonc_varid = netcdf.defVar(nc, 'lonc', 'NC_FLOAT', ...
    elem_dimid);
netcdf.putAtt(nc, lonc_varid, 'units', 'degrees_east');
netcdf.putAtt(nc, lonc_varid, 'standard_name', 'longitude');
netcdf.putAtt(nc, lonc_varid, 'long_name', 'zonal longitude');

latc_varid = netcdf.defVar(nc, 'latc', 'NC_FLOAT', ...
    elem_dimid);
netcdf.putAtt(nc, latc_varid, 'units', 'degrees_north');
netcdf.putAtt(nc, latc_varid, 'standard_name', 'latitude');
netcdf.putAtt(nc, latc_varid, 'long_name', 'zonal latitude');

zeta_varid = netcdf.defVar(nc, 'zeta', 'NC_FLOAT', ...
    [node_dimid, time_dimid]);
netcdf.putAtt(nc, zeta_varid, 'long_name', 'Water Surface Elevation');
netcdf.putAtt(nc, zeta_varid, 'units', 'meters');
netcdf.putAtt(nc, zeta_varid, 'positive', 'up');
netcdf.putAtt(nc, zeta_varid, 'standard_name', ...
    'sea_surface_height_above_geoid');
netcdf.putAtt(nc, zeta_varid, 'grid', 'Bathymetry_Mesh');
netcdf.putAtt(nc, zeta_varid, 'coordinates', 'time lat lon');
netcdf.putAtt(nc, zeta_varid, 'type', 'data');
netcdf.putAtt(nc, zeta_varid, 'location', 'node');

ua_varid = netcdf.defVar(nc, 'ua', 'NC_FLOAT', ...
    [elem_dimid, time_dimid]);
netcdf.putAtt(nc, ua_varid, 'long_name', 'Vertically Averaged x-velocity');
netcdf.putAtt(nc, ua_varid, 'units', 'meters  s-1');
netcdf.putAtt(nc, ua_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, ua_varid, 'type', 'data');

va_varid = netcdf.defVar(nc, 'va', 'NC_FLOAT', ...
    [elem_dimid, time_dimid]);
netcdf.putAtt(nc, va_varid, 'long_name', 'Vertically Averaged y-velocity');
netcdf.putAtt(nc, va_varid, 'units', 'meters  s-1');
netcdf.putAtt(nc, va_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, va_varid, 'type', 'data');

u_varid = netcdf.defVar(nc, 'u', 'NC_FLOAT', ...
    [elem_dimid, siglay_dimid, time_dimid]);
netcdf.putAtt(nc, u_varid, 'long_name', 'Eastward Water Velocity');
netcdf.putAtt(nc, u_varid, 'units', 'meters  s-1');
netcdf.putAtt(nc, u_varid, 'standard_name', 'eastward_sea_water_velocity');
netcdf.putAtt(nc, u_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, u_varid, 'coordinates', 'time siglay latc lonc');
netcdf.putAtt(nc, u_varid, 'type', 'data');
netcdf.putAtt(nc, u_varid, 'location', 'face');

v_varid = netcdf.defVar(nc, 'v', 'NC_FLOAT', ...
    [elem_dimid, siglay_dimid, time_dimid]);
netcdf.putAtt(nc, v_varid, 'long_name', 'Northward Water Velocity');
netcdf.putAtt(nc, v_varid, 'units', 'meters  s-1');
netcdf.putAtt(nc, v_varid, 'standard_name', ...
    'Northward_sea_water_velocity');
netcdf.putAtt(nc, v_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, v_varid, 'coordinates', 'time siglay latc lonc');
netcdf.putAtt(nc, v_varid, 'type', 'data');
netcdf.putAtt(nc, v_varid, 'location', 'face');

temp_varid = netcdf.defVar(nc, 'temp', 'NC_FLOAT', ...
    [node_dimid, siglay_dimid, time_dimid]);
netcdf.putAtt(nc, temp_varid, 'long_name', 'Temperature');
netcdf.putAtt(nc, temp_varid, 'standard_name', 'sea_water_temperature');
netcdf.putAtt(nc, temp_varid, 'units', 'degrees Celcius');
netcdf.putAtt(nc, temp_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, temp_varid, 'coordinates', 'time siglay lat lon');
netcdf.putAtt(nc, temp_varid, 'type', 'data');
netcdf.putAtt(nc, temp_varid, 'location', 'node');

salinity_varid = netcdf.defVar(nc, 'salinity', 'NC_FLOAT', ...
    [node_dimid, siglay_dimid, time_dimid]);
netcdf.putAtt(nc, salinity_varid, 'long_name', 'Salinity');
netcdf.putAtt(nc, salinity_varid, 'standard_name', 'sea_water_salinity');
netcdf.putAtt(nc, salinity_varid, 'units', '1e-3');
netcdf.putAtt(nc, salinity_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, salinity_varid, 'coordinates', 'time siglay lat lon');
netcdf.putAtt(nc, salinity_varid, 'type', 'data');
netcdf.putAtt(nc, salinity_varid, 'location', 'node');

hyw_varid = netcdf.defVar(nc, 'hyw', 'NC_FLOAT', ...
    [node_dimid, siglev_dimid, time_dimid]);
netcdf.putAtt(nc, hyw_varid, 'long_name', ...
    'hydro static vertical velocity');
netcdf.putAtt(nc, hyw_varid, 'units', 'meters s-1');
netcdf.putAtt(nc, hyw_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, hyw_varid, 'type', 'data');
netcdf.putAtt(nc, hyw_varid, 'coordinates', 'time siglay lat lon');

if nesttype > 2
    cweights_varid = netcdf.defVar(nc, 'weight_cell', 'NC_FLOAT', ...
    [elem_dimid, time_dimid]);
    netcdf.putAtt(nc, cweights_varid, 'long_name', ...
    'Weights for elements in relaxation zone');
    netcdf.putAtt(nc, cweights_varid, 'units', 'no units');
    netcdf.putAtt(nc, cweights_varid, 'grid', 'fvcom_grid');
    netcdf.putAtt(nc, cweights_varid, 'type', 'data');

    nweights_varid = netcdf.defVar(nc, 'weight_node', 'NC_FLOAT', ...
    [node_dimid, time_dimid]);
    netcdf.putAtt(nc, nweights_varid, 'long_name', ...
    'Weights for nodes in relaxation zone');
    netcdf.putAtt(nc, nweights_varid, 'units', 'no units');
    netcdf.putAtt(nc, nweights_varid, 'grid', 'fvcom_grid');
    netcdf.putAtt(nc, nweights_varid, 'type', 'data');
end

% end definitions
netcdf.endDef(nc);

% write time data
nStringOut = char();
[nYr, nMon, nDay, nHour, nMin, nSec] = mjulian2greg(nest.time);
for i = 1:ntimes
    if strcmp(sprintf('%02i', nSec), '60')
        % Fix some weirdness with mjulian2greg. I think this is caused by
        % rounding errors. My testing suggests this is not a problem around
        % midnight, so the number of days (and thus possibly months and
        % years) is unaffected.
        if mod(nMin(i) + 1, 60) == 0
            % Up the hour by one too
            nHour(i) = mod(nHour(i) + 1, 24);
        end
        nMin(i) = mod(nMin(i) + 1, 60);
        nSec(i) = 0;
    end
    nDate = [nYr(i), nMon(i), nDay(i), nHour(i), nMin(i), nSec(i)];
    nStringOut = [nStringOut, sprintf('%-026s', datestr(datenum(nDate), 'yyyy-mm-dd HH:MM:SS.FFF'))];
end
netcdf.putVar(nc, itime_varid, 0, numel(nest.time), floor(nest.time));
netcdf.putVar(nc, itime2_varid, 0, numel(nest.time), ...
    mod(nest.time, 1)*24*3600*1000);
netcdf.putVar(nc, Times_varid, nStringOut);
netcdf.putVar(nc, time_varid, nest.time);

% write grid information
netcdf.putVar(nc, nv_varid, nest.nv);
netcdf.putVar(nc, x_varid, nest.x);
netcdf.putVar(nc, y_varid, nest.y);
netcdf.putVar(nc, xc_varid, nest.xc);
netcdf.putVar(nc, yc_varid, nest.yc);
netcdf.putVar(nc, lon_varid, nest.lon);
netcdf.putVar(nc, lat_varid, nest.lat);
netcdf.putVar(nc, lonc_varid, nest.lonc);
netcdf.putVar(nc, latc_varid, nest.latc);

% dump data
netcdf.putVar(nc, zeta_varid, nest.zeta);
netcdf.putVar(nc, ua_varid, nest.ua);
netcdf.putVar(nc, va_varid, nest.va);
netcdf.putVar(nc, u_varid, nest.u);
netcdf.putVar(nc, v_varid, nest.v);
netcdf.putVar(nc, temp_varid, nest.temp);
netcdf.putVar(nc, salinity_varid, nest.salinity);
netcdf.putVar(nc, hyw_varid, nest.hyw);
if nesttype > 2
    netcdf.putVar(nc, cweights_varid, nest.weight_cell);
    netcdf.putVar(nc, nweights_varid, nest.weight_node);
end

% close file
netcdf.close(nc)

if ftbverbose
    fprintf('end   : %s\n', subname)
end
