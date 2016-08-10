function write_FVCOM_TT3(turbine,filename,mytitle)

% Dump tidal turbine parameters to FVCOM forcing file.
%
% function write_FVCOM_TT(turbine,filename,mytitle)
%
% DESCRIPTION:
%    Generate a NetCDF file containing tidal turbine parameters for FVCOM.
% netCDF variable names:
% turbine_sigma_layer, turbine_number, swept_area
%
% INPUT:
%   turbine.numbers = user defined array of integers specify the number of turbines in each element
%                 defined on the elements
%   turbine.sigma_layer = user defined array specifying the fractional
%   split across sigma layers that the turbine rotor swept area occupies
%   turbine.area    = rotor swept area (m^2)
%   filename  = filename to dump to
%   mytitle   = title of the case (set as global attribute)
%
% EXAMPLE USAGE:
%    write_FVCOM_TT(turbine, 'tst_file.nc', 'tst tidal turbine parameters')
%
% Author(s):
%    Rory O'Hara Murray (MSS)
%
% Revision history
%
%==============================================================================
%warning off
subname = 'write_FVCOM_TT';
global ftbverbose;
if(ftbverbose);
  fprintf('\n'); fprintf(['begin : ' subname '\n']);
end;

% check dimensions
nElems = numel(turbine.numbers);
nSiglay = size(turbine.sigma_layer,2);
if(nElems == 0)
    error('dimension of turbine_numbers is 0, something is wrong ')
end;

if isfield(turbine, 'thrust')
    write_thrust = true;
else
    write_thrust = false;
end

%------------------------------------------------------------------------------
% Dump to turbine NetCDF file
%------------------------------------------------------------------------------
if(ftbverbose);
  fprintf('Dumping to turbine parameters NetCDF file: %s\n',filename);
  fprintf('Size of turbine_numbers array: %i\n',nElems);
end;
nc = netcdf.create(filename,'clobber');

netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title',mytitle)
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history', 'File created using write_FVCOM_TT.m from the MATLAB fvcom-toolbox')

% dimensions
nele_dimid=netcdf.defDim(nc,'nele',nElems);
nsiglay_dimid=netcdf.defDim(nc,'siglay', nSiglay);

% variables and attributes
num_tt_varid=netcdf.defVar(nc,'turbine_number','NC_FLOAT',nele_dimid);
netcdf.putAtt(nc,num_tt_varid,'long_name','Number of turbines in each element');

sigma_layer_varid=netcdf.defVar(nc,'turbine_sigma_layer','NC_FLOAT',[nele_dimid nsiglay_dimid]);
netcdf.putAtt(nc,sigma_layer_varid,'long_name','the fraction of each sigma layer that the turbines are located in');

swept_area_varid=netcdf.defVar(nc,'swept_area','NC_FLOAT',nele_dimid);
netcdf.putAtt(nc,swept_area_varid,'long_name','Area swept by turbine rotor blades');

if write_thrust
    thrust_coeff_varid=netcdf.defVar(nc,'thrust_coeff','NC_FLOAT',nele_dimid);
    netcdf.putAtt(nc,thrust_coeff_varid,'long_name','Turbine thrust ceofficient');
end

% variables that could be added in the future
% blade_coeff_varid=netcdf.defVar(nc,'blade_coeff','NC_FLOAT',nele_dimid);
% netcdf.putAtt(nc,blade_coeff_varid,'long_name','Drag coefficient for the turbine blades');

% end definitions
netcdf.endDef(nc);

% write data
% netCDF variable names:
% turbine_sigma_layer, turbine_number, swept_area
netcdf.putVar(nc,num_tt_varid,turbine.numbers);
netcdf.putVar(nc,sigma_layer_varid,turbine.sigma_layer);
netcdf.putVar(nc,swept_area_varid,turbine.area);
if write_thrust
    netcdf.putVar(nc,thrust_coeff_varid,turbine.thrust);
end

% variables that could be added in the future
% netcdf.putVar(nc,blade_coeff_varid,turbine.drag);

% close file
netcdf.close(nc);

if(ftbverbose);
  fprintf(['end   : ' subname '\n'])
end;


