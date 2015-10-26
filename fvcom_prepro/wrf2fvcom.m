function wrf2fvcom(wrffile, fvfile)
% Dumps all WRF output to FVCOM compatible netCDF.
%
% wrf2fvcom(ncfile, wrffile)
%
% DESCRIPTION:
%   Take WRF netCDF output (e.g. from version 3.7.1 of WRF) and convert to
%   the format FVCOM requires for netCDF forcing files.
%
% INPUT:
%   wrffile - path to the WRF netCDF file.
%   fvfile - path to the FVCOM forcing file to create.
%
% OUTPUT:
%   FVCOM forcing netCDF file (fvfile).
%
% EXAMPLE USAGE:
%   wrf2fvcom('wrfout_d01_2000-01-01_00:00:00', 'casename_v01_wnd.nc')
%
% NOTES:
%   This currently only supports regularly gridded output (i.e. WRF native
%   grid). Support for writing out unstructured data may come in the
%   future.
%
% Author(s):
%   Dmitry Aleynik (Scottish Association for Marine Science)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2015-10-23 - First version developed with reference to Dmitry Aleynik's
%   dump_meteo_fvgrid.m function.
%
%==========================================================================

global ftbverbose;
subname = 'wrf2fvcom';
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
    fprintf('Read WRF output from %s... ', wrffile)
end

% Read the WRF output and build a struct with the relevant fields in it.
wrf_fields = {'Shortwave' ,'Net_Heat', ...
    'U10', 'V10', ...
    'Precipitation', 'Evaporation', 'SLP', ...
    'XLAT', 'XLONG', 'Times'};

for ff = 1:length(wrf_fields)
    wrf.(wrf_fields{ff}) = ncread(wrffile, wrf_fields{ff});
end

% Get the required dimensions.
finfo = ncinfo(wrffile);
dimnames = {finfo.Dimensions.Name};
for dd = 1:length(dimnames);
    switch dimnames{dd}
        case 'west_east'
            wrf.west_east = finfo.Dimensions(dd).Length;
        case 'south_north'
            wrf.south_north = finfo.Dimensions(dd).Length;
        case 'DateStrLen'
            wrf.DateStrLen = finfo.Dimensions(dd).Length;
    end
end

% Sort out times.
wrf.ntimes = length(wrf.Times(1, :));
wrf.mtime = datenum(wrf.Times');

if ftbverbose
    fprintf('done.\n')
end

grid_name = 'wrf_grid';
grid_source = 'wrf grid (structured) surface forcing';
CoordVar = 'lat lon';

fv.xlon = zeros(wrf.south_north, wrf.west_east);
fv.xlat = zeros(wrf.south_north, wrf.west_east);
fv.fmask = zeros(wrf.south_north, wrf.west_east);

[fv.xlon, fv.xlat] = deal(wrf.XLONG, wrf.XLAT);

% -------------------------------------------------------------------------
% Define the netCDF parameters
%--------------------------------------------------------------------------

if ftbverbose
    fprintf('Export WRF to FVCOM input file %s... ', fvfile)
end

% Define global attributes.
nc.type = 'FVCOM METEO FORCING FILE';
nc.title = 'WRF model forcing';
nc.history = sprintf('FILE CREATED using wrf2fvcom.m from the MATLAB fvcom-toolbox on %s', datestr(now, 31));
nc.source =  grid_source;
nc.START_DATE = datestr(wrf.mtime(1, 1), 31);
nc.END_DATE = datestr(wrf.mtime(end, 1), 31);
nc.file = fvfile;
nc.Conventions = 'CF-1.0';

% Populate the global attributes.
ncid = netcdf.create(fvfile, 'clobber');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'type', nc.type);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', nc.title);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'history', nc.history);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'source', nc.source);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'START_DATE', nc.START_DATE);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'END_DATE', nc.END_DATE);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'file', nc.file);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'Conventions',nc.Conventions);

sn_idim = netcdf.defDim(ncid, 'south_north', wrf.south_north);
we_idim = netcdf.defDim(ncid, 'west_east', wrf.west_east);
DateStrLen_idim = netcdf.defDim(ncid, 'DateStrLen', wrf.DateStrLen);
time_idim = netcdf.defDim(ncid, 'time', 0);

% Variables
varid_Times = netcdf.defVar(ncid, 'Times', 'NC_CHAR', [DateStrLen_idim, time_idim]);
netcdf.putAtt(ncid, varid_Times, 'description', 'UTC');

varid_XLAT = netcdf.defVar(ncid, 'XLAT', 'NC_FLOAT', [we_idim, sn_idim]);
varid_XLONG = netcdf.defVar(ncid, 'XLONG', 'NC_FLOAT', [we_idim, sn_idim]);
varid_Shortwave = netcdf.defVar(ncid, 'Shortwave', 'NC_FLOAT', [we_idim, sn_idim, time_idim]);
varid_Net_Heat = netcdf.defVar(ncid, 'Net_Heat', 'NC_FLOAT', [we_idim, sn_idim, time_idim]);
varid_U10 = netcdf.defVar(ncid, 'U10', 'NC_FLOAT', [we_idim, sn_idim, time_idim]);
varid_V10 = netcdf.defVar(ncid, 'V10', 'NC_FLOAT', [we_idim, sn_idim, time_idim]);
varid_Precipitation = netcdf.defVar(ncid, 'Precipitation', 'NC_FLOAT', [we_idim, sn_idim, time_idim]);
varid_Evaporation = netcdf.defVar(ncid, 'Evaporation', 'NC_FLOAT', [we_idim, sn_idim, time_idim]);
varid_SLP = netcdf.defVar(ncid, 'air_pressure', 'NC_FLOAT',[we_idim, sn_idim, time_idim]);

netcdf.putAtt(ncid, varid_XLAT, 'long_name', 'latitude');
netcdf.putAtt(ncid, varid_XLAT, 'description', 'LATITUDE, SOUTH IS NEGATIVE');
netcdf.putAtt(ncid, varid_XLAT, 'units', 'degrees_north');
netcdf.putAtt(ncid, varid_XLAT, 'type', 'data');
netcdf.putAtt(ncid, varid_XLONG, 'long_name', 'longitude');
netcdf.putAtt(ncid, varid_XLONG, 'description', 'LONGITUDE, WEST IS NEGATIVE');
netcdf.putAtt(ncid, varid_XLONG, 'units', 'degrees_east');
netcdf.putAtt(ncid, varid_XLONG, 'type', 'data');

netcdf.putAtt(ncid, varid_Shortwave, 'long_name', 'Shortwave, upward is negative');
netcdf.putAtt(ncid, varid_Shortwave, 'units', 'W m-2');
netcdf.putAtt(ncid, varid_Shortwave, 'grid', grid_name);
netcdf.putAtt(ncid, varid_Shortwave, 'coordinates', CoordVar);
netcdf.putAtt(ncid, varid_Shortwave, 'type', 'data');

netcdf.putAtt(ncid, varid_Net_Heat, 'long_name', 'Sum of shortwave, longwave, sensible and latent heat fluxes, ocean lose heat is negative');
netcdf.putAtt(ncid, varid_Net_Heat, 'units', 'W m-2');
netcdf.putAtt(ncid, varid_Net_Heat, 'grid', grid_name);
netcdf.putAtt(ncid, varid_Net_Heat, 'coordinates', CoordVar);
netcdf.putAtt(ncid, varid_Net_Heat, 'type', 'data');

netcdf.putAtt(ncid, varid_U10, 'long_name', 'Eastward Wind Velocity');
netcdf.putAtt(ncid, varid_U10, 'description', 'U at 10 M');
netcdf.putAtt(ncid, varid_U10, 'units', 'm s-1'  );
netcdf.putAtt(ncid, varid_U10, 'grid', grid_name);
netcdf.putAtt(ncid, varid_U10, 'type', 'data');

netcdf.putAtt(ncid, varid_V10, 'long_name', 'Northward Wind Velocity');
netcdf.putAtt(ncid, varid_V10, 'description', 'V at 10 M');
netcdf.putAtt(ncid, varid_V10, 'units', 'm s-1'  );
netcdf.putAtt(ncid, varid_V10, 'grid', grid_name);
netcdf.putAtt(ncid, varid_V10, 'type', 'data');

netcdf.putAtt(ncid, varid_Precipitation, 'long_name', 'Precipitation');
netcdf.putAtt(ncid, varid_Precipitation, 'description', 'Precipitation, ocean lose water is negative');
netcdf.putAtt(ncid, varid_Precipitation, 'units', 'm s-1'  );
netcdf.putAtt(ncid, varid_Precipitation, 'grid', grid_name);
netcdf.putAtt(ncid, varid_Precipitation, 'coordinates', CoordVar);
netcdf.putAtt(ncid, varid_Precipitation, 'type', 'data');

netcdf.putAtt(ncid, varid_Evaporation, 'long_name', 'Evaporation');
netcdf.putAtt(ncid, varid_Evaporation, 'description', 'Evaporation, ocean lose water is negative');
netcdf.putAtt(ncid, varid_Evaporation, 'units', 'm s-1'  );
netcdf.putAtt(ncid, varid_Evaporation, 'grid', grid_name);
netcdf.putAtt(ncid, varid_Evaporation, 'coordinates', CoordVar);
netcdf.putAtt(ncid, varid_Evaporation, 'type', 'data');

netcdf.putAtt(ncid, varid_SLP, 'long_name', 'Air Pressure');
netcdf.putAtt(ncid, varid_SLP, 'description', 'Sea surface airpressure');
netcdf.putAtt(ncid, varid_SLP, 'units', 'Pa'  );
netcdf.putAtt(ncid, varid_SLP, 'grid', grid_name);
netcdf.putAtt(ncid, varid_SLP, 'coordinates', CoordVar);
netcdf.putAtt(ncid, varid_SLP, 'type', 'data');

netcdf.endDef(ncid);

% -------------------------------------------------------------------------
% Dump the data
% -------------------------------------------------------------------------

for i = 1:wrf.ntimes
    cc = datestr(wrf.mtime(i), 31);
    netcdf.putVar(ncid, varid_Times, [0, i-1], [19, 1], cc);
end

netcdf.putVar(ncid, varid_Shortwave, [0, 0, 0], size(wrf.Shortwave), wrf.Shortwave); % shortwave W/m^2
netcdf.putVar(ncid, varid_Net_Heat, [0, 0, 0], size(wrf.Net_Heat), wrf.Net_Heat); % net W/m^2
netcdf.putVar(ncid, varid_U10, [0, 0, 0], size(wrf.U10), wrf.U10); % m/s
netcdf.putVar(ncid, varid_V10, [0, 0, 0], size(wrf.V10), wrf.V10); % m/s
netcdf.putVar(ncid, varid_Precipitation, [0, 0, 0], size(wrf.Precipitation), wrf.Precipitation); % m/s
netcdf.putVar(ncid, varid_Evaporation, [0, 0, 0], size(wrf.Evaporation), wrf.Evaporation); % m/s
netcdf.putVar(ncid, varid_SLP, [0, 0, 0], size(wrf.SLP), wrf.SLP * 100); % mbar -> Pa

netcdf.close(ncid);

if ftbverbose
    fprintf('done.\n')
end

if ftbverbose
    fprintf('end   : %s \n', subname)
end