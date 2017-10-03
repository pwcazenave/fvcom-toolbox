function write_FVCOM_sediment(sediment, filename)
% Write spatially-varying sediment threshold data.
%
% function write_FVCOM_sediment(sediment, filename)
%
% DESCRIPTION:
%    Generate a netCDF file containing spatially variable sediment propoerties
%    for FVCOM
%
% INPUT
%   sediment  = struct with nodal data in the following mandatory fields:
%               - tau_ce
%               - tau_cd (this appears to be unused but required at launch)
%               - erate
%               and the following optional field(s):
%               - <sediment_class>_bedfrac
%                If the optional field(s) is included, the namelist option
%                SEDIMENT_PARAMETER_TYPE = "non-uniform" can be used.
%   filename  = filename to which to write the outputs.
%
% OUTPUT:
%    netCDF file: filename of the output netCDF file to generate.
%
% EXAMPLE USAGE
%    write_FVCOM_sediment(sediment, 'tst_sediment.nc')
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2017-03-27 New function to create a netCDF for the SEDIMENT_PARAMETER_FILE
%    section in the model namelist.
%
%==========================================================================

global ftbverbose
[~, subname] = fileparts(mfilename('fullpath'));
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

%--------------------------------------------------------------------------
% Check inputs.
%--------------------------------------------------------------------------
assert(exist('sediment', 'var') == 1, 'incorrect usage of %s, must provide sediment struct', subname)
assert(exist('filename', 'var') == 1, 'incorrect usage of %s, must provide ''filename''', subname)

% Check dimensions
nVerts = numel(sediment.tau_cd);
assert(nVerts ~= 0, 'dimension of sediment data is 0, something is wrong')

%--------------------------------------------------------------------------
% Dump to netCDF file
%--------------------------------------------------------------------------
if ftbverbose
    fprintf('Writing sediment data to netCDF file: %s\n', filename);
    fprintf('nodes dimension: %d\n', nVerts);
end

nc = netcdf.create(filename, 'clobber');

netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'), 'history', ...
    sprintf('File created using %s from the MATLAB fvcom-toolbox', subname))


% define dimensions
node_dimid = netcdf.defDim(nc, 'node', nVerts);

% define variables and attributes
erate_varid = netcdf.defVar(nc, 'erate', 'NC_DOUBLE', node_dimid);
netcdf.putAtt(nc, erate_varid, 'long_name', 'bed erodibility constant');
netcdf.putAtt(nc, erate_varid, 'units', '-');

tau_ce_varid = netcdf.defVar(nc, 'tau_ce', 'NC_DOUBLE', node_dimid);
netcdf.putAtt(nc, tau_ce_varid, 'long_name', 'critical shear stress for erosion');
netcdf.putAtt(nc, tau_ce_varid, 'units', 'N/m^2');

tau_cd_varid = netcdf.defVar(nc, 'tau_cd', 'NC_DOUBLE', node_dimid);
netcdf.putAtt(nc, tau_cd_varid, 'long_name', 'critical shear stress for deposition');
netcdf.putAtt(nc, tau_cd_varid, 'units', 'N/m^2');

extra_fields = setdiff(fieldnames(sediment), {'tau_ce', 'tau_cd', 'erate'});
for f = 1:length(extra_fields)
    sediment_name = extra_fields{f};
    if ftbverbose
        fprintf('Defining extra variable: %s\n', sediment_name)
    end
    eval(sprintf('%s_varid = netcdf.defVar(nc, ''%s'', ''NC_DOUBLE'', node_dimid);\n', sediment_name, sediment_name))
    eval(sprintf('netcdf.putAtt(nc, %s_varid, ''long_name'', ''Fraction of %s in surface layer'');\n', sediment_name, sediment_name))
    eval(sprintf('netcdf.putAtt(nc, %s_varid, ''units'', ''-'')\n', sediment_name))
    eval(sprintf('netcdf.putAtt(nc, %s_varid, ''grid'', ''fvcom_grid'')\n', sediment_name))
end

% end definitions
netcdf.endDef(nc);

% dump data
netcdf.putVar(nc, erate_varid, sediment.erate);
netcdf.putVar(nc, tau_ce_varid, sediment.tau_ce);
netcdf.putVar(nc, tau_cd_varid, sediment.tau_cd);

% Add the extra data.
for f = 1:length(extra_fields)
    sediment_name = extra_fields{f};
    if ftbverbose
        fprintf('Adding extra variable: %s\n', sediment_name)
    end
    eval(sprintf('netcdf.putVar(nc, %s_varid, sediment.(sediment_name));', sediment_name))
end

% close netCDF
netcdf.close(nc)

if ftbverbose
    fprintf('end   : %s\n', subname)
end
