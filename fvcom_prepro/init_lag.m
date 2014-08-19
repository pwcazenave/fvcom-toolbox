function init_lag(Mobj, timerange, type, file)
% Write online Lagrangian initial conditions file.
%
% init_lag(Mobj, type, file)
%
% DESCRIPTION:
%   Create the Lagrangian initial positions netCDF file.
%
% INPUT:
%   Mobj - MATLAB mesh object containing fields:
%       have_bath - boolean flag indicating whether the bathymetry is
%       loaded.
%       native_coord - string of the native coordinate type.
%       h - water depth on nodes.
%       xc, yc or lonc, lat - element centre coordinates (cartesian or
%       spherical depending on native_coord).
%   timerange - start and end date array (in Modified Julian Days).
%   type - supply an array of element centre indices.
%   file - path to the output netCDF file (string).
%
% OUTPUT:
%   FVCOM forcing netCDF file(s)
%
% EXAMPLE USAGE:
%   init = '/path/to/output/casename_lag.nc';
%   init_lag(Mobj, 1, init);
%
% TODO:
%   Add pathlength, group and mark as arguments to the function (currently
%   set to 0, 1 and 0 respectively).
%
% Author(s):
%   Geoff Cowles (University of Massachusetts Dartmouth)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2014-08-19 First version based on Geoff's original example_init_lag.m
%   file.
%
%==========================================================================

subname = 'init_lag';

global ftbverbose;
if ftbverbose
	fprintf('\nbegin : %s \n', subname)
end

tbeg = timerange(1);
tend = timerange(2);

if type == 1 % initialize at all elements
    xp   = Mobj.xc;
    yp   = Mobj.yc;
    nLag = Mobj.nElems;
elseif type == 2 % initialize along a line of interest
    nLag = 10;
    p1 = [1.188363e6,194497];
    p2 = [1.188548e6,194996];
    xp = p1(1):(p2(1)-p1(1))/(nLag-1):p2(1);
    yp = p1(2):(p2(2)-p1(2))/(nLag-1):p2(2);
elseif isvector(type)
end

% Dump the initial particle position file.
nc = netcdf.create(file, 'clobber');
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'references', 'http://fvcom.smast.umassd.edu, http://pml.ac.uk')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'source', 'init_lag.m')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'info', 'debugging')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'history', 'File generated using init_lag.m from the MATLAB fvcom-toolbox')

% Dimensions
nlag_dimid = netcdf.defDim(nc, 'nparticles', nLag);

% Particle vars
x_varid = netcdf.defVar(nc, 'x', 'NC_FLOAT', nlag_dimid);
netcdf.putAtt(nc, x_varid, 'long_name', 'particle x position');
netcdf.putAtt(nc, x_varid, 'units', 'm');

y_varid = netcdf.defVar(nc, 'particle y position', 'NC_FLOAT', nlag_dimid);
netcdf.putAtt(nc, y_varid, 'long_name', 'particle y position');
netcdf.putAtt(nc, y_varid, 'units', 'm');

z_varid = netcdf.defVar(nc, 'particle z position', 'NC_FLOAT', nlag_dimid);
netcdf.putAtt(nc, z_varid, 'long_name', 'particle z position');
netcdf.putAtt(nc, z_varid, 'units', 'm');

pathlength_varid = netcdf.defVar(nc, 'pathlength', 'NC_FLOAT', nlag_dimid);
netcdf.putAtt(nc, pathlength_varid, 'long_name', 'particle integrated path length');
netcdf.putAtt(nc, pathlength_varid, 'units', 'm');

tbeg_varid = netcdf.defVar(nc, 'tbeg', 'NC_FLOAT', nlag_dimid);
netcdf.putAtt(nc, tbeg_varid, 'long_name', 'particle release time');
netcdf.putAtt(nc, tbeg_varid, 'units', 'days since 1858-11-17 00:00:00');
netcdf.putAtt(nc, tbeg_varid, 'format', 'modified julian day (MJD)');
netcdf.putAtt(nc, tbeg_varid, 'time_zone', 'UTC');

tend_varid = netcdf.defVar(nc, 'tbeg', 'NC_FLOAT', nlag_dimid);
netcdf.putAtt(nc, tend_varid, 'long_name', 'particle freeze time');
netcdf.putAtt(nc, tend_varid, 'units', 'days since 1858-11-17 00:00:00');
netcdf.putAtt(nc, tend_varid, 'format', 'modified julian day (MJD)');
netcdf.putAtt(nc, tend_varid, 'time_zone', 'UTC');

group_varid = netdf.defVar(nc, 'group', 'NC_INT', nlag_dimid);
netcdf.putAtt(nc, group_varid, 'long_name', 'particle group');
netcdf.putAtt(nc, group_varid, 'units', '-');

mark_varid = netdf.defVar(nc, 'mark', 'NC_INT', nlag_dimid);
netcdf.putAtt(nc, mark_varid, 'long_name', 'particle mark');
netcdf.putAtt(nc, mark_varid, 'units', '-');

% End definitions
netcdf.endDef(nc);

% Dump vars
netcdf.putVar(nc, x_varid, xp);
netcdf.putVar(nc, y_varid, yp);
netcdf.putVar(nc, z_varid, zp);
netcdf.putVar(nc, tbeg_varid, tbeg);
netcdf.putVar(nc, tend_varid, tend);
netcdf.putVar(nc, group_varid, 1);
netcdf.putVar(nc, mark_varid, 0);
netcdf.putVar(nc, pathlength_varid, 0);

% Close netCDF file
netcdf.close(nc);

if ftbverbose; fprintf('end   : %s\n', subname); end
