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
%
% Revision history
%   
%==============================================================================
warning off
global ftbverbose;
if(ftbverbose);
  subname = 'write_FVCOM_bedflag';
  fprintf('\n'); fprintf(['begin : ' subname '\n']);
end;

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
if(~exist('bedflag'))
	error('incorrect usage of gen_bedflag_file, must provide bedflag field')
end;
if(~exist('filename'))
	error('incorrect usage of gen_bedflag_file, must provide filename')
end;
if(~exist('title'))
	error('incorrect usage of gen_bedflag_file, must provide title field')
end;

% check dimensions
nVerts = prod(size(bedflag));
if(nVerts == 0)
	error('dimension of bedflag is 0, something is wrong ')
end;

%------------------------------------------------------------------------------
% Dump to bedflag NetCDF file
%------------------------------------------------------------------------------
if(ftbverbose);
fprintf('Dumping to bedflag NetCDF file: \n',filename);
fprintf('Size of bedflag array: \n',nVerts);
end;
nc = netcdf(filename,'clobber');
nc.title = mytitle;
nc('node') = prod(size(bedflag));
nc{'bedflag'}  = ncfloat('node');
nc{'bedflag'}.long_name = 'bed deposition flag';
nc{'bedflag'}.units = '-';
nc{'bedflag'}(1:nVerts) = bedflag(1:nVerts);
ierr = close(nc);



if(ftbverbose); fprintf(['end   : ' subname '\n']); end;


