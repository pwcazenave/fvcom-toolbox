function write_FVCOM_bedflag(bedflag,filename,mytitle) 

% Dump spatially-variable flag (bedflag) to FVCOM forcing file
%
% function write_FVCOM_bedflag(bedflag,filename,mytitle)
%
% DESCRIPTION:
%    Generate a NetCDF file containing spatially variable bedflag for FVCOM 
%
% INPUT 
%   bedflag   = user defined bed flag (=0, no erosion/bedosition, =1, erosion/bedosition) 
%               on the nodes
%   filename  = filename to dump to
%   mytitle   = title of the case (set as global attribute) 
%
% OUTPUT:
%    NetCDF file: filename
%
% EXAMPLE USAGE
%    write_FVCOM_bedflag(bedflag, 'tst_bedflag.nc', 'no bedosition/erosion in Skagit river')
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2016-02-18 Updated the code to use the MATLAB netCDF routines.
%   
%==============================================================================

global ftbverbose
subname = 'write_FVCOM_bedflag';
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
if(~exist('bedflag'))
	error('incorrect usage of gen_bedflag_file, must provide bedflag field')
end
if(~exist('filename'))
	error('incorrect usage of gen_bedflag_file, must provide filename')
end
if(~exist('title'))
	error('incorrect usage of gen_bedflag_file, must provide title field')
end

% check dimensions
nVerts = numel(bedflag);
if(nVerts == 0)
	error('dimension of bedflag is 0, something is wrong ')
end;

%------------------------------------------------------------------------------
% Dump to bedflag NetCDF file
%------------------------------------------------------------------------------
if ftbverbose
    fprintf('Dumping to bedflag NetCDF file: %s\m', filename);
    fprintf('Size of bedflag array: %d\n', nVerts);
end

nc = netcdf.create(filename, 'clobber');

% define dimensions
node_dimid=netcdf.defDim(nc, 'node', numel(bedflag));

% define variables and attributes
node_varid=netcdf.defVar(nc, 'bedflag', 'NC_INT', node_dimid);
netcdf.putAtt(nc,node_varid, 'long_name', 'bed deposition flag');
netcdf.putAtt(nc,node_varid, 'units', '-');

% end definitions
netcdf.endDef(nc);

% dump data
netcdf.putVar(nc, node_varid, bedflag);

% close netCDF
netcdf.close(nc)

if ftbverbose
    fprintf('end   : %s', subname)
end
