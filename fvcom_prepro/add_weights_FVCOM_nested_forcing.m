function add_weights_FVCOM_nested_forcing(conf, ncfile2read)
% Creates an FVCOM nesting file.
%
% function add_weights_FVCOM_nested_forcing(conf, ncfile)
%
% DESCRIPTION:
%   Modifies in place the nested file to add weights for the number of levels 
%   suitable for types 3. 
%
%
% INPUT:
%   conf        = with either new_weights field or conf.nest.levels
%   ncfile2read = full path to the nesting file to be created.
%
% OUTPUT:
%   FVCOM nesting file.
%
% EXAMPLE USAGE:
%   conf.new_weight_cell = weights see manual for explanation [0-1]. It
%           only requires nlevel values, not spatial or temporal dimension is
%           expected.
%   conf.new_weight_node = weights see manual for explanation [0-1]. It
%           only requires nlevel values, not spatial or temporal dimension is
%           expected.
%   conf.nest.levels = [1]; The value reflect the inner most level to keep.
%           1 will keep the most external boundary, 4 will keep 1-4 levels.
%
%   modify_FVCOM_nested_forcing(conf, '/tmp/fvcom_nested.nc', 1)
%
% Author(s):
%   Ricardo Torres (Plymouth Marine Laboratory)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2017-01-17 First version based on Pierre's write_FVCOM_nested_forcing.m
%   script. Compress the time series data to save space. Requires
%   netCDF4 in FVCOM.
%
%==========================================================================

% The following variables are read from the file:
%
% lon, lat:     Grid node positions             [node]
% lonc, latc:   Grid element positions          [nele]
% h:            Grid node depths                [node]
% hc:           Grid element depth              [nele]
% nv:           Triangulation table             [nele, 3]
% zeta:         Sea surface elevation           [node, time]
% ua:           Vertically averaged x velocity  [node, time]
% va:           Vertically averaged y velocity  [nele, time]
% u:            Eastward Water Velocity         [nele, siglay, time]
% v:            Northward Water Velocity        [nele, siglay, time]
% temp:         Temperature                     [node, siglay, time]
% salinity:     Salinity                        [node, siglay, time]
% hyw:          Hydrostatic vertical velocity   [node, siglev, time]
% weight_cell:  Weighting for elements          [nele]
% weight_node:  Weighting for nodes             [node]
% Itime:        Days since 1858-11-17 00:00:00  [time]
% Itime2:       msec since 00:00:00             [time]
%
% We include these optional ones for humans:
% time:         Modified Julian Day             [time]
% Times:        Gregorian dates                 [time, datestrlen]

[~, subname] = fileparts(mfilename('fullpath'));

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

% Can't use CLOBBER and NETCDF4 at the same time (the bitwise or didn't
% work). Fall back to a horrible delete and then create instead.
if exist(ncfile2read, 'file')
    nc2read = netcdf.open(ncfile2read,'WRITE');
    [PATHSTR,NAME,EXT] = fileparts(ncfile2read);
    if ftbverbose
        fprintf('\nProcessing file : %s\n', ncfile2read)
    end
else
    sprintf('I am stoping. The nesting file %s doesn''t exist',ncfile2read);
    error(['Aborting the script ',subname])
end
% List of data we need.
required = {'weight_cell', 'weight_node'};
for f = required
    nest.(f{1})=[];
end
% Read all variables from the existing file
if ~isfield (conf,'nest')
    error('Missing conf.nest, aborting... ')
end
if ~isfield (conf.nest,'levels')
    error('Missing conf.nest.levels, aborting... ')
end
% define dimensions

timeid = netcdf.inqDimID(nc2read, 'time');
siglayid = netcdf.inqDimID(nc2read, 'siglay');
siglevid = netcdf.inqDimID(nc2read, 'siglev');
nodeid = netcdf.inqDimID(nc2read,'node');
neleid = netcdf.inqDimID(nc2read,'nele');
       
[~,nsiglay] = netcdf.inqDim(nc2read, siglayid);
[~,nsiglev] = netcdf.inqDim(nc2read, siglevid);
[~,ntimes] = netcdf.inqDim(nc2read, timeid);
[~,nodes] = netcdf.inqDim(nc2read, nodeid);
[~,elems] = netcdf.inqDim(nc2read, neleid);

% Replicate weights in time

nest.weight_cell = repmat([conf.new_weight_cell{:}]',1,ntimes);
nest.weight_node = repmat([conf.new_weight_node{:}]',1,ntimes)

[elem_nest] =size(nest.weight_cell,1);

if elem_nest ~= elems
    error('We have different elements in FVCOM nesting output file and our calculations')
end
try
    netcdf.reDef(nc2read);

    cweights_varid = netcdf.defVar(nc2read, 'weight_cell', 'NC_FLOAT', ...
        [neleid, timeid]);
    netcdf.putAtt(nc2read, cweights_varid, 'long_name', ...
        'Weights for elements in relaxation zone');
    netcdf.putAtt(nc2read, cweights_varid, 'units', 'no units');
    netcdf.putAtt(nc2read, cweights_varid, 'grid', 'fvcom_grid');
    netcdf.putAtt(nc2read, cweights_varid, 'type', 'data');
    
    nweights_varid = netcdf.defVar(nc2read, 'weight_node', 'NC_FLOAT', ...
        [nodeid, timeid]);
    netcdf.putAtt(nc2read, nweights_varid, 'long_name', ...
        'Weights for nodes in relaxation zone');
    netcdf.putAtt(nc2read, nweights_varid, 'units', 'no units');
    netcdf.putAtt(nc2read, nweights_varid, 'grid', 'fvcom_grid');
    netcdf.putAtt(nc2read, nweights_varid, 'type', 'data');

% end definitions
netcdf.endDef(nc2read);


    netcdf.putVar(nc2read, cweights_varid, nest.weight_cell);
    netcdf.putVar(nc2read, nweights_varid, nest.weight_node);
catch e
    fprintf(e.message)
    error('Adding variable %s failed - does the variable already exist?', 'weight_cell')
end


% close file
netcdf.close(nc2read)

if ftbverbose
    fprintf('end   : %s\n', subname)
end
