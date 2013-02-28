function write_FVCOM_meanflow(Mobj, ncfile, data)
% Export mean flow forcing at the open boundary to NetCDF.
%
% function write_FVCOM_meanflow(Mobj, ncfile, data)
%
% DESCRIPTION:
%    Setup an FVCOM hydrographic open boundary mean flow forcing file.
%
% INPUT:
%   Mobj    - MATLAB mesh object (with fields mf_time, siglay, siglev,
%               nObcElements and read_obc_elements).
%   ncfile  - Output NetCDF file name.
%   data    - 2D array of mean flow along the open boundary sized
%               (nobcelems, time).
%
% OUTPUT:
%    FVCOM mean flow values along the FVCOM open boundary in a NETCDF file
%    named ncfile.
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% TODO: 
%    Add support for multiple open boundaries (currently hard limit of
%    one).
%
% Revision history
%    2013-02-20 - First version.
%    2013-02-28 - Change output of velocities to be at the boundary element
%    centres rather than the boundary nodes.
%
%==========================================================================

subname = 'write_FVCOM_meanflow';

global ftbverbose
if ftbverbose
    fprintf('\n'); fprintf(['begin : ' subname '\n']);
end

% open new NetCDF file
nc = netcdf.create(ncfile, 'clobber');

% define dimensions
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'), 'type', 'FVCOM MEANFLOW TIME SERIES FILE')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'), 'title', 'FVCOM MEANFLOW TIME SERIES data for open boundary')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'), 'history', 'File created using the write_FVCOM_meanflow.m from the MATLAB fvcom-toolbox')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'), 'filename', ncfile)
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'), 'Conventions', 'CF-1.0')

% define dimensions
nmfcell_dimid = netcdf.defDim(nc, 'nmfcell', Mobj.nObcElements);
time_dimid = netcdf.defDim(nc, 'time', netcdf.getConstant('NC_UNLIMITED'));
siglay_dimid = netcdf.defDim(nc, 'siglay', size(Mobj.siglay, 2));
siglev_dimid = netcdf.defDim(nc, 'siglev', size(Mobj.siglev, 2));

time_varid = netcdf.defVar(nc, 'time', 'NC_FLOAT', time_dimid);
netcdf.putAtt(nc, time_varid, 'long_name', 'time');
netcdf.putAtt(nc, time_varid, 'units', 'days since 1858-11-17 00:00:00');
netcdf.putAtt(nc, time_varid, 'format', 'modified julian day (MJD)');
netcdf.putAtt(nc, time_varid, 'time_zone', 'UTC');

itime_varid = netcdf.defVar(nc, 'Itime', 'NC_INT', time_dimid);
netcdf.putAtt(nc, itime_varid, 'units', 'days since 1858-11-17 00:00:00');
netcdf.putAtt(nc, itime_varid, 'format', 'modified julian day (MJD)');
netcdf.putAtt(nc, itime_varid, 'time_zone', 'UTC');

itime2_varid = netcdf.defVar(nc, 'Itime2', 'NC_INT', time_dimid);
netcdf.putAtt(nc, itime2_varid, 'units', 'msec since 00:00:00');
netcdf.putAtt(nc, itime2_varid, 'time_zone', 'UTC');

nmfcell_varid = netcdf.defVar(nc, 'I_MFCELL_GL', 'NC_INT', nmfcell_dimid);
netcdf.putAtt(nc, nmfcell_varid, 'long_name', 'Open Boundary Cell Number');
netcdf.putAtt(nc, nmfcell_varid, 'grid', 'obc_grid');
netcdf.putAtt(nc, nmfcell_varid, 'type', 'data');

mfdist_varid = netcdf.defVar(nc, 'MFDIST', 'NC_FLOAT', [siglay_dimid, nmfcell_dimid]);
netcdf.putAtt(nc, mfdist_varid, 'long_name', 'Mean Flow Flux Vertical Distribution');
netcdf.putAtt(nc, mfdist_varid, 'units', 'no units');
netcdf.putAtt(nc, mfdist_varid, 'grid', 'obc_grid');
netcdf.putAtt(nc, mfdist_varid, 'type', 'data');

dmfqdis_varid = netcdf.defVar(nc, 'DMFQDIS', 'NC_FLOAT', [nmfcell_dimid, time_dimid]);
netcdf.putAtt(nc, dmfqdis_varid, 'long_name', 'open boundary mean flow flux');
netcdf.putAtt(nc, dmfqdis_varid, 'units', 'm^3/s');
netcdf.putAtt(nc, dmfqdis_varid, 'grid', 'obc_grid');
netcdf.putAtt(nc, dmfqdis_varid, 'type', 'data');

% end definitions
netcdf.endDef(nc);

% write data
netcdf.putVar(nc, time_varid, 0, numel(Mobj.mf_times), Mobj.mf_times);
netcdf.putVar(nc, itime_varid, floor(Mobj.mf_times));
netcdf.putVar(nc, itime2_varid, 0, numel(Mobj.mf_times), mod(Mobj.mf_times, 1) * 24 * 3600 * 1000);
netcdf.putVar(nc, nmfcell_varid, Mobj.read_obc_elements{1});
% MFDIST is calculated here as the diff of the sigma levels. This should
% work for uniform and gaussian etc. distributions so long as Mobj.siglev
% has the right values in it.
netcdf.putVar(nc, mfdist_varid, repmat(abs(diff(Mobj.siglev)), [numel([Mobj.read_obc_elements{:}]), 1])');
netcdf.putVar(nc, dmfqdis_varid, [0, 0], [numel(Mobj.read_obc_elements{1}), numel(Mobj.mf_times)], data);

netcdf.close(nc);

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end

