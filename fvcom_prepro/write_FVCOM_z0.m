function write_FVCOM_z0(z0,filename,mytitle,cbcmin)
% Dump spatially-variable or uniform bottom roughness (z0) to FVCOM forcing
% file.
%
% function write_FVCOM_z0(z0, filename, mytitle)
%
% DESCRIPTION:
%    Generate a NetCDF file containing spatially variable z0 for FVCOM
%
% INPUT
%   z0        = user defined roughness field (m)
%               roughness is defined on the elements
%               expect values between 3 10^-3 (gravel) and .2 10^-3 (i.e. 0.0002) for mud
%   filename  = filename to dump to
%   mytitle   = title of the case (set as global attribute)
%   cbcmin    = minimum value of CBC (optional, defaults to 0.0018 if
%               omitted).
%
% OUTPUT:
%    netCDF file called `filename'
%
% EXAMPLE USAGE
%    write_FVCOM_z0(z0field, 'tst_z0.nc', 'z0 tst domain')
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2012-06-15 Added support for native MATLAB NetCDF routines. Requires
%    MATLAB 2010a or higher.
%    2016-08-09 Added new variable (cbcmin) to support FVCOM version 4.
%    Also tidied up the code a bit.
%
%==============================================================================

[~, subname] = fileparts(mfilename('fullpath'));
global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s\n', subname);
end

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
if ~exist('z0', 'var')
    error('incorrect usage of write_FVCOM_z0, must provide z0 field')
end
if ~exist('filename', 'var')
    error('incorrect usage of write_FVCOM_z0, must provide filename')
end
if ~exist('mytitle', 'var')
    error('incorrect usage of write_FVCOM_z0, must provide title field')
end

% check dimensions
nElems = numel(z0);
if nElems == 0
    error('Number of elements in z0 is 0.')
end

% If we haven't been given a value of cbc min, set it to the example from
% Jianzhong Ge.
if nargin == 3
    cbcmin = repmat(0.0018, nElems, 1);
end

%------------------------------------------------------------------------------
% Dump to variables to the netCDF file
%------------------------------------------------------------------------------
if ftbverbose
  fprintf('Dumping to z0 NetCDF file: %s\n', filename);
  fprintf('Size of z0 array: %i\n', nElems);
end
nc = netcdf.create(filename, 'clobber');

netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'title', mytitle)
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'history', ...
    sprintf('File created with %s from the MATLAB fvcom-toolbox', subname))

% dimensions
nele_dimid = netcdf.defDim(nc, 'nele', nElems);

% variables and attributes
z0b_varid = netcdf.defVar(nc, 'z0b', 'NC_FLOAT', nele_dimid);
netcdf.putAtt(nc, z0b_varid, 'long_name', 'bottom roughness');
netcdf.putAtt(nc, z0b_varid, 'units', 'm');
netcdf.putAtt(nc, z0b_varid, 'type', 'data');

cbcmin_varid=netcdf.defVar(nc, 'cbcmin', 'NC_FLOAT', nele_dimid);
netcdf.putAtt(nc, cbcmin_varid, 'long_name', 'bottom roughness minimum');
netcdf.putAtt(nc, cbcmin_varid, 'units', 'None');
netcdf.putAtt(nc, cbcmin_varid, 'type', 'data');

% end definitions
netcdf.endDef(nc);

% write data
netcdf.putVar(nc, z0b_varid, z0);
netcdf.putVar(nc, cbcmin_varid, cbcmin);

% close file
netcdf.close(nc);

if ftbverbose
  fprintf('end   : %s\n', subname)
end
