function write_FVCOM_groundwater(Mobj, grndwtr_file, varargin)
% Write an FVCOM groundwater time series forcing file.
%
% write_FVCOM_groundwater(Mobj, grndwtr_file, varargin)
%
% DESCRIPTION:
%   Write an FVCOM groundwater time series forcing file.
%
% INPUT:
%   Mobj = Matlab mesh object with fields:
%       obc_nodes - array of boundary node IDs.
%       surfaceElevation - array of surface elevation values (shaped [space,
%       time]).
%   grndwtr_file = name of netCDF file.
%
%   Optional keyword-argument pairs. These control the time variables. This
%   script defaults to writing 'Times' only.
%   FVCOM needs only one of:
%       1. Times: character string of times
%       2. Itime and Itime2: integer days and milliseconds since midnight
%       3. time: float days.
%   FVCOM checks for these in the order above and this script defaults to
%   writing Times only. Adjust the keyword-argument pairs to your liking:
%
%   'strtime' = set to true to output the 'Times' variable
%   'inttime' = set to true to output the 'Itime' and 'Itime2' variables
%   'floattime' = set to true to output the 'time' variable
%
% OUTPUT:
%   grndwtr_file = a netCDF FVCOM surface elevations tide forcing file.
%
% EXAMPLE USAGE
%   With default settings:
%       write_FVCOM_groundwater(Mobj, '/tmp/grndwtr.nc)
%   Enable the 'time' variable in the netCDF.
%       write_FVCOM_groundwater(Mobj, '/tmp/grndwtr.nc, ...
%           'floattime', true)
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2015-11-18 First version.
%
%==========================================================================

global ftbverbose

subname = 'write_FVCOM_groundwater';
if ftbverbose
    fprintf('\nbegin : %s \n', subname);
end

% Default to string times as FVCOM looks for these first.
strtime = true;
inttime = false;
floattime = false;
for vv = 1:2:length(varargin)
    switch varargin{vv}
        case 'strtime'
            strtime = true;
        case 'inttime'
            inttime = true;
        case 'floattime'
            floattime = true;
    end
end

%%
%--------------------------------------------------------------------------
% Dump the file
%--------------------------------------------------------------------------

nc = netcdf.create(grndwtr_file, 'clobber');

% Define global attributes
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'type', 'FVCOM GROUNDWATER FORCING FILE')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'source', 'FVCOM grid (unstructured) surface forcing')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'history', sprintf('File created with %s from the MATLAB fvcom-toolbox', subname))

% Dfine dimensions
node_dimid = netcdf.defDim(nc, 'node', Mobj.nVerts);
% Although the nele dimension isn't used, FVCOM checks for its existence
% and FATAL_ERRORs if it's missing, so we need to add it here.
nele_dimid = netcdf.defDim(nc, 'nele', Mobj.nElems); %#ok<NASGU>
time_dimid = netcdf.defDim(nc, 'time', netcdf.getConstant('NC_UNLIMITED'));
date_str_len_dimid = netcdf.defDim(nc, 'DateStrLen', 26);

% Define variables and attributes
gflux_varid = netcdf.defVar(nc, 'groundwater_flux', 'NC_FLOAT', [node_dimid, time_dimid]);
netcdf.putAtt(nc, gflux_varid, 'long_name', 'groundwater volume flux');
netcdf.putAtt(nc, gflux_varid, 'units', 'm3 s-1');
netcdf.putAtt(nc, gflux_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, gflux_varid, 'type', 'data');

gtemp_varid = netcdf.defVar(nc, 'groundwater_temp', 'NC_FLOAT', [node_dimid, time_dimid]);
netcdf.putAtt(nc, gtemp_varid, 'long_name', 'groundwater inflow temperature');
netcdf.putAtt(nc, gtemp_varid, 'units', 'degrees_C');
netcdf.putAtt(nc, gtemp_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, gtemp_varid, 'type', 'data');

gsalt_varid = netcdf.defVar(nc, 'groundwater_salt', 'NC_FLOAT', [node_dimid, time_dimid]);
netcdf.putAtt(nc, gsalt_varid, 'long_name', 'groundwater inflow salinity');
netcdf.putAtt(nc, gsalt_varid, 'units', '1e-3');
netcdf.putAtt(nc, gsalt_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, gsalt_varid, 'type', 'data');

iint_varid = netcdf.defVar(nc, 'iint', 'NC_INT', time_dimid);
netcdf.putAtt(nc, iint_varid, 'long_name', 'internal mode iteration number');

if floattime
    time_varid = netcdf.defVar(nc, 'time', 'NC_FLOAT', time_dimid);
    netcdf.putAtt(nc, time_varid, 'long_name', 'time');
    netcdf.putAtt(nc, time_varid, 'units', 'days since 1858-11-17 00:00:00');
    netcdf.putAtt(nc, time_varid, 'format', 'modified julian day (MJD)');
    netcdf.putAtt(nc, time_varid, 'time_zone', 'UTC');
end

if inttime
    itime_varid = netcdf.defVar(nc, 'Itime', 'NC_INT', time_dimid);
    netcdf.putAtt(nc, itime_varid, 'units', 'days since 1858-11-17 00:00:00');
    netcdf.putAtt(nc, itime_varid, 'format', 'modified julian day (MJD)');
    netcdf.putAtt(nc, itime_varid, 'time_zone', 'UTC');

    itime2_varid = netcdf.defVar(nc, 'Itime2', 'NC_INT', time_dimid);
    netcdf.putAtt(nc, itime2_varid, 'units', 'msec since 00:00:00');
    netcdf.putAtt(nc, itime2_varid, 'time_zone', 'UTC');
end

if strtime
    Times_varid = netcdf.defVar(nc, 'Times', 'NC_CHAR', [date_str_len_dimid, time_dimid]);
    netcdf.putAtt(nc, Times_varid, 'time_zone', 'UTC');
end

% End definitions
netcdf.endDef(nc);

% Write data
nTimes = length(Mobj.groundwater.times);
netcdf.putVar(nc, iint_varid, 0, nTimes, 1:nTimes);
if floattime
    netcdf.putVar(nc, time_varid, 0, nTimes, Mobj.groundwater.times);
end
if inttime
    netcdf.putVar(nc, itime_varid, floor(Mobj.groundwater.times));
    netcdf.putVar(nc, itime2_varid, 0, nTimes, round(mod(Mobj.groundwater.times, 1) * 24 * 60 * 60 * 1000));
end
if strtime
    nStringOut = char();
    [nYr, nMon, nDay, nHour, nMin, nSec] = mjulian2greg(Mobj.groundwater.times);
    for i = 1:nTimes
        nDate = [nYr(i),  nMon(i),  nDay(i),  nHour(i),  nMin(i),  nSec(i)];
        nStringOut = [nStringOut, sprintf('%04i/%02i/%02i %02i:%02i:%09.6f',  nDate)];
    end
    netcdf.putVar(nc, Times_varid, nStringOut);
end
netcdf.putVar(nc, gflux_varid, Mobj.groundwater.flux);
netcdf.putVar(nc, gtemp_varid, Mobj.groundwater.temp);
netcdf.putVar(nc, gsalt_varid, Mobj.groundwater.salt);

% Close file
netcdf.close(nc);

if ftbverbose
    fprintf('end   : %s \n',  subname)
end
