function add_var_FVCOM_river(RiverFile,VarName,VarLongName,VarUnits,VarData)
% add time dependent scalar variable to a Riverfile
%
% function add_var_FVCOM_river(RiverFile,VarName,VarLongName,VarUnits,VarData)
%
% DESCRIPTION:
%    Write an additional scalar variable (e.g. sediment, DO) to a netCDF
%    River file.  Note that the concentration of the scalar variable
%    is set the same at all river points in the file so it is assumed
%    that even if the file contains multiple nodes, they all refer to the same
%    river.
%
% INPUT
%    RiverFile:   FVCOM 3.x netCDF river forcing file
%    VarName:     Variable name (will be the name of the array in the netCDF file)
%    VarLongName: Variable attribute "long_name"
%    VarUnits:    Variable attribute "units"
%    VarData:     1-D time series of variable data of exact same dimensions as
%                 the river flux
%
% OUTPUT:
%    Modified FVCOM RiverFile
%
% EXAMPLE USAGE
%    add_var_FVCOM_river('tst_riv.nc','medium_sand','medium sand','kg-m^-3',sand_ts)
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2016-04-01 Updated the code to use the native MATLAB netCDF routines.
%
%==============================================================================

%warning off

[~, subname] = fileparts(mfilename('fullpath'));
global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

if ftbverbose
    fprintf('adding variable %s to file %s\n', VarName, RiverFile)
end

%------------------------------------------------------------------------------
% Open River Netcdf and read dimensions
%------------------------------------------------------------------------------
if exist(RiverFile, 'file') ~= 2
    error('file: %s does not exist', RiverFile)
end

% read dimensions
flux = ncread(RiverFile, 'river_flux');
[nTimes, nRivnodes] = size(flux);

% make sure time dimension matches FVCOM river dims
tempTimes = prod(size(VarData));
if nTimes ~= tempTimes
    fprintf('# of time frames in file %s is %d\n', RiverFile, tempTimes)
    fprintf('# of time frames in VarData is %d\n', nTimes)
    error('You have chosen the wrong vocation')
end


%------------------------------------------------------------------------------
% Write variable definition and data and close file
%------------------------------------------------------------------------------

% set field
river_field = repmat(VarData, nRivnodes, 1);

%------------------------------------------------------------------------------
% dump to netcdf file
%------------------------------------------------------------------------------

% open boundary forcing
nc = netcdf.open(RiverFile, 'NC_WRITE');

% define dimensions
time_dimid = netcdf.inqDimID(nc, 'time');
rivers_dimid = netcdf.inqDimID(nc, 'rivers');

% add the new variable if it doesn't already exist.
try
    netcdf.reDef(nc);
    varid = netcdf.defVar(nc, VarName, 'NC_FLOAT', [rivers_dimid, time_dimid]);
    netcdf.putAtt(nc, varid, 'long_name', VarLongName);
    netcdf.putAtt(nc, varid, 'units', VarUnits);

    % end definitions
    netcdf.endDef(nc);

    % dump dynamic data
    netcdf.putVar(nc, varid, river_field);
catch e
    fprintf(e.message)
    error('Adding variable %s failed - does the variable already exist?', VarName)
end

netcdf.close(nc);

if ftbverbose
    fprintf('end   : %s\n', subname)
end

