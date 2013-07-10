function write_FVCOM_z0(z0,filename,mytitle)

% Dump spatially-variable or uniform bottom roughness (z0) to FVCOM forcing
% file.
%
% function write_FVCOM_z0(z0,filename,mytitle)
%
% DESCRIPTION:
%    Generate a NetCDF file containing spatially variable z0 for FVCOM
%
% INPUT
%   z0        = user defined roughness field (m)
%               roughness is defined on the elements
%   filename  = filename to dump to
%   mytitle   = title of the case (set as global attribute)
%
% OUTPUT:
%    NetCDF file: filename
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
%
%==============================================================================
warning off
subname = 'write_FVCOM_z0';
global ftbverbose;
if(ftbverbose);
  fprintf('\n'); fprintf(['begin : ' subname '\n']);
end;

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
if(~exist('z0','var'))
	error('incorrect usage of write_FVCOM_z0, must provide z0 field')
end;
if(~exist('filename','var'))
	error('incorrect usage of write_FVCOM_z0, must provide filename')
end;
if(~exist('mytitle','var'))
	error('incorrect usage of write_FVCOM_z0, must provide title field')
end;

% check dimensions
nElems = numel(z0);
if(nElems == 0)
	error('dimension of z0 is 0, something is wrong ')
end;

%------------------------------------------------------------------------------
% Dump to z0 NetCDF file
%------------------------------------------------------------------------------
if(ftbverbose);
  fprintf('Dumping to z0 NetCDF file: %s\n',filename);
  fprintf('Size of z0 array: %i\n',nElems);
end;
nc = netcdf.create(filename,'clobber');

netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title',mytitle)

% dimensions
nele_dimid=netcdf.defDim(nc,'nele',nElems);

% variables and attributes
z0b_varid=netcdf.defVar(nc,'z0b','NC_FLOAT',nele_dimid);
netcdf.putAtt(nc,z0b_varid,'long_name','bottom roughness');
netcdf.putAtt(nc,z0b_varid,'units','m');

% end definitions
netcdf.endDef(nc);

% write data
netcdf.putVar(nc,z0b_varid,z0);

% close file
netcdf.close(nc);

if(ftbverbose);
  fprintf(['end   : ' subname '\n'])
end;


