function modify_FVCOM_nested_forcing(conf, ncfile2read, nesttype)
% Creates an FVCOM nesting file.
%
% function modify_FVCOM_nested_forcing(conf, ncfile, nesttype)
%
% DESCRIPTION:
%   Modifies nested file to change weights or number of levels and writes a
%   new Nested file suitable for types different from 3. you can only do
%   one or the other. Changing weights and reducing the number of levels is
%   not permitted in the current implementation. You can always cut the
%   number of levels and then call it again to change the values of the
%   weights.
%
%
% INPUT:
%   conf        = with either new_weights field or conf.nest.levels
%   ncfile2read = full path to the nesting file to be created.
%   nesttype    = [optional] nesting type (defaults to 1 = direct nesting).
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
if nargin == 2
    nesttype = 1;
elseif nargin < 2 || nargin > 3
    error(['Incorrect input arguments. Supply netCDF file path, ', ...
        'conf struct and optionally the nesting type (1, 2 or 3).'])
end

% Can't use CLOBBER and NETCDF4 at the same time (the bitwise or didn't
% work). Fall back to a horrible delete and then create instead.
if exist(ncfile2read, 'file')
    nc2read = netcdf.open(ncfile2read,'WRITE');
    [PATHSTR,NAME,EXT] = fileparts(ncfile2read);
    % new file to hold the modified nesting data
    ncfile = fullfile(PATHSTR,[NAME,'modified',EXT]);
else
    sprintf('I am stoping. The nesting file %s doesn''t exist',ncfile2read);
    error(['Aborting the script ',subname])
end
% List of data we need.
required = {'time', 'x', 'y', 'lon', 'lat', 'xc', 'yc', 'lonc', 'latc', ...
    'nv','zeta', 'h', 'h_center', 'u', 'v', 'ua', 'va', 'temp', 'salinity', 'hyw', ...
    'weight_cell', 'weight_node', 'siglay', 'siglay_center', 'siglev', 'siglev_center'};
for f = required
    nest.(f{1})=[];
end
% If we only need to change the weights we don't need to rewrite the entire
% file
if isfield(conf,'new_weight_cell') && isfield(conf,'new_weight_node')
    varid = netcdf.inqVarID(nc2read,'weight_cell');
    nest.weight_cell = netcdf.getVar(nc2read,varid,'double');
    weight_cell = unique(conf.new_weight_cell); % these are sorted
    Nweight_cell = unique(nest.weight_cell); % these are sorted
    if length(weight_cell) ~= length(Nweight_cell)
        error('new_weight_cell doesn''t have the same number of UNIQUE levels as the nesting file')
    end
    for rr=1:length(Nweight_cell)
        nest.weight_cell(nest.weight_cell==Nweight_cell(rr)) = weight_cell(rr);
    end
    netcdf.putVar(nc2read, varid, nest.weight_cell);
    % now do the Node's weights
    varid = netcdf.inqVarID(nc2read,'weight_node');
    nest.weight_node = netcdf.getVar(nc2read,varid,'double');
    weight_node = unique(conf.new_weight_node); % these are sorted
    Nweight_node = unique(nest.weight_node); % these are sorted
    if length(weight_node) ~= length(Nweight_node)
        error('new_weight_node doesn''t have the same number of UNIQUE levels as the nesting file')
    end
    for rr=1:length(Nweight_node)
        nest.weight_node(nest.weight_node==Nweight_node(rr)) = weight_node(rr);
    end
    netcdf.putVar(nc2read, varid, nest.weight_node);
    sprintf('The weights have been changed in the original nesting file  %s ',ncfile2read);
    disp('Returning to the main script...')

    return

end
% Read all variables from the existing file
if ~isfield (conf,'nest')
            error('Missing conf.nest, aborting... ')
end
if ~isfield (conf.nest,'levels')
            error('Missing conf.nest.levels, aborting... ')
end

for f = required
    try
        varid = netcdf.inqVarID(nc2read,f{1});
    catch
        error( 'Missing %s variable from nesting file', f{1});
    end

    nest.(f{1}) = netcdf.getVar(nc2read,varid,'double');
end
[~, nsiglay, ~] = size(nest.u);
nsiglev = nsiglay + 1;
[~, ~] = size(nest.zeta);

% Split nodes into the different levels
    Nweight_node = unique(nest.weight_node); % these are sorted increasing
    Nweight_cell = unique(nest.weight_cell); % these are sorted increasing
    nlevs = unique(conf.nest.levels);
% Decide how many levels to keep
warning(['Chopping levels from nesting file. Assuming unique weight values for',...
    'each nesting level. If there are repeated values, this function will not work'])
    levels2keepN = Nweight_node(end:-1:end-nlevs+1); % these are now decreasing in order
    levels2keepC = Nweight_cell(end:-1:end-nlevs+1); % these are now decreasing in order
% remove nodes and cells from each variable
nodeid = netcdf.inqDimID(nc2read,'node');
neleid = netcdf.inqDimID(nc2read,'nele');
dropidxnode = nest.weight_node(:,1)<min(levels2keepN);
dropidxele = nest.weight_cell(:,1)<min(levels2keepC);
for f = required
        varid = netcdf.inqVarID(nc2read,f{1});
        [~,~,dimids,~] = netcdf.inqVar(nc2read,varid);
% check the variable has a node or element dimension
        if any(dimids==nodeid)
            switch ndims (nest.(f{1}))
                case 1
            nest.(f{1})(dropidxnode)=[];
                case 2
            nest.(f{1})(dropidxnode,:)=[];

                case 3
            nest.(f{1})(dropidxnode,:,:)=[];
            end
        elseif any(dimids==neleid)
            switch ndims (nest.(f{1}))
                case 1
            nest.(f{1})(dropidxele)=[];
                case 2
            nest.(f{1})(dropidxele,:)=[];

                case 3
            nest.(f{1})(dropidxele,:,:)=[];
            end

        end
end

% Can't use CLOBBER and NETCDF4 at the same time (the bitwise or didn't
% work). Fall back to a horrible delete and then create instead.
if exist(ncfile, 'file')
    delete(ncfile)
end
[elems, nsiglay, ntimes] = size(nest.u);
nsiglev = nsiglay + 1;
[nodes, ~] = size(nest.zeta);

% output the new nesting file
nc = netcdf.create(ncfile, 'NETCDF4');

% define global attributes
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'type', ...
    'FVCOM nesting TIME SERIES FILE')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'title', ...
    sprintf('FVCOM nesting TYPE %d TIME SERIES data for open boundary', ...
    nesttype))
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'history', ...
    sprintf('File created using %s from the MATLAB fvcom-toolbox', subname))
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'filename', ncfile)
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'Conventions', 'CF-1.0')

% define dimensions
elem_dimid = netcdf.defDim(nc, 'nele', elems);
node_dimid = netcdf.defDim(nc, 'node', nodes);
three_dimid = netcdf.defDim(nc, 'three', 3);
time_dimid = netcdf.defDim(nc, 'time', netcdf.getConstant('NC_UNLIMITED'));
siglay_dimid = netcdf.defDim(nc, 'siglay', nsiglay);
siglev_dimid = netcdf.defDim(nc, 'siglev', nsiglev);
datestrlen_dimid = netcdf.defDim(nc, 'DateStrLen', 26);

% define variables
time_varid = netcdf.defVar(nc, 'time', 'NC_FLOAT', time_dimid);
netcdf.putAtt(nc, time_varid, 'long_name', 'time');
netcdf.putAtt(nc, time_varid, 'units', 'days since 1858-11-17 00:00:00');
netcdf.putAtt(nc, time_varid, 'format', 'modified julian day (MJD)');
netcdf.putAtt(nc, time_varid, 'time_zone', 'UTC');

itime_varid = netcdf.defVar(nc, 'Itime', 'NC_INT', ...
    time_dimid);
netcdf.putAtt(nc, itime_varid, 'units', 'days since 1858-11-17 00:00:00');
netcdf.putAtt(nc, itime_varid, 'format', 'modified julian day (MJD)');
netcdf.putAtt(nc, itime_varid, 'time_zone', 'UTC');

itime2_varid = netcdf.defVar(nc, 'Itime2', 'NC_INT', ...
    time_dimid);
netcdf.putAtt(nc, itime2_varid, 'units', 'msec since 00:00:00');
netcdf.putAtt(nc, itime2_varid, 'time_zone', 'UTC');

Times_varid = netcdf.defVar(nc, 'Times' ,'NC_CHAR', ...
    [datestrlen_dimid, time_dimid]);
netcdf.putAtt(nc, Times_varid, 'time_zone', 'UTC');

x_varid = netcdf.defVar(nc, 'x', 'NC_FLOAT', ...
    node_dimid);
netcdf.putAtt(nc, x_varid, 'units', 'meters');
netcdf.putAtt(nc, x_varid, 'long_name', 'nodal x-coordinate');

y_varid = netcdf.defVar(nc, 'y', 'NC_FLOAT', ...
    node_dimid);
netcdf.putAtt(nc, y_varid, 'units', 'meters');
netcdf.putAtt(nc, y_varid, 'long_name', 'nodal y-coordinate');

xc_varid = netcdf.defVar(nc, 'xc', 'NC_FLOAT', ...
    elem_dimid);
netcdf.putAtt(nc, xc_varid, 'units', 'meters');
netcdf.putAtt(nc, xc_varid, 'long_name', 'zonal x-coordinate');

yc_varid = netcdf.defVar(nc, 'yc', 'NC_FLOAT', ...
    elem_dimid);
netcdf.putAtt(nc, yc_varid, 'units', 'meters');
netcdf.putAtt(nc, yc_varid, 'long_name', 'zonal y-coordinate');

lon_varid = netcdf.defVar(nc, 'lon', 'NC_FLOAT', ...
    node_dimid);
netcdf.putAtt(nc, lon_varid, 'units', 'degrees_east');
netcdf.putAtt(nc, lon_varid, 'standard_name', 'longitude');
netcdf.putAtt(nc, lon_varid, 'long_name', 'nodal longitude');

lat_varid = netcdf.defVar(nc, 'lat', 'NC_FLOAT', ...
    node_dimid);
netcdf.putAtt(nc, lat_varid, 'units', 'degrees_north');
netcdf.putAtt(nc, lat_varid, 'standard_name', 'latitude');
netcdf.putAtt(nc, lat_varid, 'long_name', 'nodal latitude');

lonc_varid = netcdf.defVar(nc, 'lonc', 'NC_FLOAT', ...
    elem_dimid);
netcdf.putAtt(nc, lonc_varid, 'units', 'degrees_east');
netcdf.putAtt(nc, lonc_varid, 'standard_name', 'longitude');
netcdf.putAtt(nc, lonc_varid, 'long_name', 'zonal longitude');

latc_varid = netcdf.defVar(nc, 'latc', 'NC_FLOAT', ...
    elem_dimid);
netcdf.putAtt(nc, latc_varid, 'units', 'degrees_north');
netcdf.putAtt(nc, latc_varid, 'standard_name', 'latitude');
netcdf.putAtt(nc, latc_varid, 'long_name', 'zonal latitude');

nv_varid = netcdf.defVar(nc, 'nv', 'NC_INT', ...
    [elem_dimid, three_dimid]);
netcdf.putAtt(nc, xc_varid, 'units', 'no units');
netcdf.putAtt(nc, xc_varid, 'long_name', 'elements nodes indices');

zeta_varid = netcdf.defVar(nc, 'zeta', 'NC_FLOAT', ...
    [node_dimid, time_dimid]);
netcdf.putAtt(nc, zeta_varid, 'long_name', 'Water Surface Elevation');
netcdf.putAtt(nc, zeta_varid, 'units', 'meters');
netcdf.putAtt(nc, zeta_varid, 'positive', 'up');
netcdf.putAtt(nc, zeta_varid, 'standard_name', ...
    'sea_surface_height_above_geoid');
netcdf.putAtt(nc, zeta_varid, 'grid', 'Bathymetry_Mesh');
netcdf.putAtt(nc, zeta_varid, 'coordinates', 'time lat lon');
netcdf.putAtt(nc, zeta_varid, 'type', 'data');
netcdf.putAtt(nc, zeta_varid, 'location', 'node');

ua_varid = netcdf.defVar(nc, 'ua', 'NC_FLOAT', ...
    [elem_dimid, time_dimid]);
netcdf.putAtt(nc, ua_varid, 'long_name', 'Vertically Averaged x-velocity');
netcdf.putAtt(nc, ua_varid, 'units', 'meters  s-1');
netcdf.putAtt(nc, ua_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, ua_varid, 'type', 'data');

va_varid = netcdf.defVar(nc, 'va', 'NC_FLOAT', ...
    [elem_dimid, time_dimid]);
netcdf.putAtt(nc, va_varid, 'long_name', 'Vertically Averaged y-velocity');
netcdf.putAtt(nc, va_varid, 'units', 'meters  s-1');
netcdf.putAtt(nc, va_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, va_varid, 'type', 'data');

u_varid = netcdf.defVar(nc, 'u', 'NC_FLOAT', ...
    [elem_dimid, siglay_dimid, time_dimid]);
netcdf.putAtt(nc, u_varid, 'long_name', 'Eastward Water Velocity');
netcdf.putAtt(nc, u_varid, 'units', 'meters  s-1');
netcdf.putAtt(nc, u_varid, 'standard_name', 'eastward_sea_water_velocity');
netcdf.putAtt(nc, u_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, u_varid, 'coordinates', 'time siglay latc lonc');
netcdf.putAtt(nc, u_varid, 'type', 'data');
netcdf.putAtt(nc, u_varid, 'location', 'face');

v_varid = netcdf.defVar(nc, 'v', 'NC_FLOAT', ...
    [elem_dimid, siglay_dimid, time_dimid]);
netcdf.putAtt(nc, v_varid, 'long_name', 'Northward Water Velocity');
netcdf.putAtt(nc, v_varid, 'units', 'meters  s-1');
netcdf.putAtt(nc, v_varid, 'standard_name', ...
    'Northward_sea_water_velocity');
netcdf.putAtt(nc, v_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, v_varid, 'coordinates', 'time siglay latc lonc');
netcdf.putAtt(nc, v_varid, 'type', 'data');
netcdf.putAtt(nc, v_varid, 'location', 'face');

temp_varid = netcdf.defVar(nc, 'temp', 'NC_FLOAT', ...
    [node_dimid, siglay_dimid, time_dimid]);
netcdf.putAtt(nc, temp_varid, 'long_name', 'Temperature');
netcdf.putAtt(nc, temp_varid, 'standard_name', 'sea_water_temperature');
netcdf.putAtt(nc, temp_varid, 'units', 'degrees Celcius');
netcdf.putAtt(nc, temp_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, temp_varid, 'coordinates', 'time siglay lat lon');
netcdf.putAtt(nc, temp_varid, 'type', 'data');
netcdf.putAtt(nc, temp_varid, 'location', 'node');

salinity_varid = netcdf.defVar(nc, 'salinity', 'NC_FLOAT', ...
    [node_dimid, siglay_dimid, time_dimid]);
netcdf.putAtt(nc, salinity_varid, 'long_name', 'Salinity');
netcdf.putAtt(nc, salinity_varid, 'standard_name', 'sea_water_salinity');
netcdf.putAtt(nc, salinity_varid, 'units', '1e-3');
netcdf.putAtt(nc, salinity_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, salinity_varid, 'coordinates', 'time siglay lat lon');
netcdf.putAtt(nc, salinity_varid, 'type', 'data');
netcdf.putAtt(nc, salinity_varid, 'location', 'node');

hyw_varid = netcdf.defVar(nc, 'hyw', 'NC_FLOAT', ...
    [node_dimid, siglev_dimid, time_dimid]);
netcdf.putAtt(nc, hyw_varid, 'long_name', ...
    'hydro static vertical velocity');
netcdf.putAtt(nc, hyw_varid, 'units', 'meters s-1');
netcdf.putAtt(nc, hyw_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, hyw_varid, 'type', 'data');
netcdf.putAtt(nc, hyw_varid, 'coordinates', 'time siglay lat lon');

siglay_varid = netcdf.defVar(nc, 'siglay', 'NC_FLOAT', ...
    [node_dimid, siglay_dimid]);
netcdf.putAtt(nc, siglay_varid, 'long_name', 'Sigma Layers');
netcdf.putAtt(nc, siglay_varid, 'standard_name', 'ocean_sigma/general_coordinate');
netcdf.putAtt(nc, siglay_varid, 'positive', 'up');
netcdf.putAtt(nc, siglay_varid, 'valid_min', '-1.f');
netcdf.putAtt(nc, siglay_varid, 'valid_max', '0.f');
netcdf.putAtt(nc, siglay_varid, 'formula_terms', 'sigma: siglay eta: zeta depth: h');

siglayc_varid = netcdf.defVar(nc, 'siglay_center', 'NC_FLOAT', ...
    [elem_dimid, siglay_dimid]);
netcdf.putAtt(nc, siglayc_varid, 'long_name', 'Sigma Layers');
netcdf.putAtt(nc, siglayc_varid, 'standard_name', 'ocean_sigma/general_coordinate');
netcdf.putAtt(nc, siglayc_varid, 'positive', 'up');
netcdf.putAtt(nc, siglayc_varid, 'valid_min', '-1.f');
netcdf.putAtt(nc, siglayc_varid, 'valid_max', '0.f');
netcdf.putAtt(nc, siglayc_varid, 'formula_terms', 'sigma: siglay_center eta: zeta_center depth: h_center');

siglev_varid = netcdf.defVar(nc, 'siglev', 'NC_FLOAT', ...
    [node_dimid, siglev_dimid]);
netcdf.putAtt(nc, siglev_varid, 'long_name', 'Sigma Levels');
netcdf.putAtt(nc, siglev_varid, 'standard_name', 'ocean_sigma/general_coordinate');
netcdf.putAtt(nc, siglev_varid, 'positive', 'up');
netcdf.putAtt(nc, siglev_varid, 'valid_min', '-1.f');
netcdf.putAtt(nc, siglev_varid, 'valid_max', '0.f');
netcdf.putAtt(nc, siglev_varid, 'formula_terms', 'sigma:siglev eta: zeta depth: h');

siglevc_varid = netcdf.defVar(nc, 'siglev_center', 'NC_FLOAT', ...
    [elem_dimid, siglev_dimid]);
netcdf.putAtt(nc, siglevc_varid, 'long_name', 'Sigma Layers');
netcdf.putAtt(nc, siglevc_varid, 'standard_name', 'ocean_sigma/general_coordinate');
netcdf.putAtt(nc, siglevc_varid, 'positive', 'up');
netcdf.putAtt(nc, siglevc_varid, 'valid_min', '-1.f');
netcdf.putAtt(nc, siglevc_varid, 'valid_max', '0.f');
netcdf.putAtt(nc, siglevc_varid, 'formula_terms', 'sigma: siglev_center eta: zeta_center depth: h_center');

h_varid = netcdf.defVar(nc, 'h', 'NC_FLOAT', ...
    node_dimid);
netcdf.putAtt(nc, h_varid, 'long_name', 'Bathymetry');
netcdf.putAtt(nc, h_varid, 'standard_name', 'sea_floor_depth_below_geoid');
netcdf.putAtt(nc, h_varid, 'units', 'm');
netcdf.putAtt(nc, h_varid, 'positive', 'down');
netcdf.putAtt(nc, h_varid, 'grid', 'Bathymetry_mesh');
netcdf.putAtt(nc, h_varid, 'coordinates', 'x y');
netcdf.putAtt(nc, h_varid, 'type', 'data');

hc_varid = netcdf.defVar(nc, 'h_center', 'NC_FLOAT', ...
    elem_dimid);
netcdf.putAtt(nc, hc_varid, 'long_name', 'Bathymetry');
netcdf.putAtt(nc, hc_varid, 'standard_name', 'sea_floor_depth_below_geoid');
netcdf.putAtt(nc, hc_varid, 'units', 'm');
netcdf.putAtt(nc, hc_varid, 'positive', 'down');
netcdf.putAtt(nc, hc_varid, 'grid', 'grid1 grid3');
netcdf.putAtt(nc, hc_varid, 'coordinates', 'latc lonc');
netcdf.putAtt(nc, hc_varid, 'grid_location', 'center');

if nesttype > 2
    cweights_varid = netcdf.defVar(nc, 'weight_cell', 'NC_FLOAT', ...
    [elem_dimid, time_dimid]);
    netcdf.putAtt(nc, cweights_varid, 'long_name', ...
    'Weights for elements in relaxation zone');
    netcdf.putAtt(nc, cweights_varid, 'units', 'no units');
    netcdf.putAtt(nc, cweights_varid, 'grid', 'fvcom_grid');
    netcdf.putAtt(nc, cweights_varid, 'type', 'data');

    nweights_varid = netcdf.defVar(nc, 'weight_node', 'NC_FLOAT', ...
    [node_dimid, time_dimid]);
    netcdf.putAtt(nc, nweights_varid, 'long_name', ...
    'Weights for nodes in relaxation zone');
    netcdf.putAtt(nc, nweights_varid, 'units', 'no units');
    netcdf.putAtt(nc, nweights_varid, 'grid', 'fvcom_grid');
    netcdf.putAtt(nc, nweights_varid, 'type', 'data');
end

% enable compression on the big variables.
netcdf.defVarDeflate(nc, zeta_varid, true, true, 7);
netcdf.defVarDeflate(nc, u_varid, true, true, 7);
netcdf.defVarDeflate(nc, v_varid, true, true, 7);
netcdf.defVarDeflate(nc, temp_varid, true, true, 7);
netcdf.defVarDeflate(nc, salinity_varid, true, true, 7);
netcdf.defVarDeflate(nc, hyw_varid, true, true, 7);

% end definitions
netcdf.endDef(nc);

% write time data
nStringOut = char();
[nYr, nMon, nDay, nHour, nMin, nSec] = mjulian2greg(nest.time);
for i = 1:ntimes
    if strcmp(sprintf('%02i', nSec), '60')
        % Fix some weirdness with mjulian2greg. I think this is caused by
        % rounding errors. My testing suggests this is not a problem around
        % midnight, so the number of days (and thus possibly months and
        % years) is unaffected.
        if mod(nMin(i) + 1, 60) == 0
            % Up the hour by one too
            nHour(i) = mod(nHour(i) + 1, 24);
        end
        nMin(i) = mod(nMin(i) + 1, 60);
        nSec(i) = 0;
    end
    nDate = [nYr(i), nMon(i), nDay(i), nHour(i), nMin(i), nSec(i)];
    nStringOut = [nStringOut, sprintf('%-026s', datestr(datenum(nDate), 'yyyy-mm-dd HH:MM:SS.FFF'))];
end
netcdf.putVar(nc, itime_varid, 0, numel(nest.time), floor(nest.time));
netcdf.putVar(nc, itime2_varid, 0, numel(nest.time), ...
    mod(nest.time, 1)*24*3600*1000);
netcdf.putVar(nc, Times_varid, nStringOut);
netcdf.putVar(nc, time_varid, nest.time);

% write grid information
netcdf.putVar(nc, nv_varid, nest.nv);
netcdf.putVar(nc, x_varid, nest.x);
netcdf.putVar(nc, y_varid, nest.y);
netcdf.putVar(nc, xc_varid, nest.xc);
netcdf.putVar(nc, yc_varid, nest.yc);
netcdf.putVar(nc, lon_varid, nest.lon);
netcdf.putVar(nc, lat_varid, nest.lat);
netcdf.putVar(nc, lonc_varid, nest.lonc);
netcdf.putVar(nc, latc_varid, nest.latc);

% dump data
netcdf.putVar(nc, zeta_varid, nest.zeta);
netcdf.putVar(nc, ua_varid, nest.ua);
netcdf.putVar(nc, va_varid, nest.va);
netcdf.putVar(nc, u_varid, nest.u);
netcdf.putVar(nc, v_varid, nest.v);
netcdf.putVar(nc, temp_varid, nest.temp);
netcdf.putVar(nc, salinity_varid, nest.salinity);
netcdf.putVar(nc, hyw_varid, nest.hyw);
netcdf.putVar(nc, siglay_varid, nest.siglay);
netcdf.putVar(nc, siglayc_varid, nest.siglay_center);
netcdf.putVar(nc, siglev_varid, nest.siglev);
netcdf.putVar(nc, siglevc_varid, nest.siglev_center);
netcdf.putVar(nc, h_varid, nest.h);
netcdf.putVar(nc, hc_varid, nest.h_center);
if nesttype > 2
    netcdf.putVar(nc, cweights_varid, nest.weight_cell);
    netcdf.putVar(nc, nweights_varid, nest.weight_node);
end

% close file
netcdf.close(nc)

if ftbverbose
    fprintf('end   : %s\n', subname)
end
