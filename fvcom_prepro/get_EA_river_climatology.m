function Mobj = get_EA_river_climatology(Mobj, ea, dist_thresh)
% Read river temperature climatologies from the Environment Agency river
% temperature data. If no data are found within the threshold specified, a
% mean climatology from all data points is provided instead.
%
% function Mobj = get_EA_river(Mogj, ea)
%
% DESCRIPTION:
%   Load all the river data from the river climatology netCDF file and find
%   the most relevant one to the river nodes in Mobj.rivers.positions.
%
% INPUT:
%   Mobj        : MATLAB mesh structure which must contain:
%                   - Mobj.river_nodes - river node IDs.
%                   - Mobj.lon, Mobj.lat - unstructured grid node
%                   positions.
%   ea          : Full path to the river climatology netCDF file.
%   dist_thresh : distance threshold beyond which a river temperature
%                 climatology data point is considered too far to be valid
%                 (in degrees).
%
% OUTPUT:
%   Mobj        : MATLAB structure with a new Mobj.river_temp field which
%                 contains the climatology for the river nodes.
%
% EXAMPLE USAGE:
%   Mobj = get_EA_river_climatology(Mobj, '/path/to/netcdf.nc', 0.05)
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Reivision history
%   2013-11-05 First version.

subname = 'get_EA_river_climatology';

global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

nc = netcdf.open(ea, 'NOWRITE');

% Load the position (lon/lat), time and climatology variables only. Not
% really bothered about the other variables.
varid = netcdf.inqVarID(nc, 'climatology');
climatology = netcdf.getVar(nc, varid, 'single');
varid = netcdf.inqVarID(nc, 'time');
time = netcdf.getVar(nc, varid, 'single');
varid = netcdf.inqVarID(nc, 'lon');
lon = netcdf.getVar(nc, varid, 'single');
varid = netcdf.inqVarID(nc, 'lat');
lat = netcdf.getVar(nc, varid, 'single');

netcdf.close(nc)

% Now find the nearest nodes to the river node positions.
nr = length(Mobj.river_nodes);

clim = nan(length(time), nr);

for r = 1:nr
    dist = sqrt((lon - Mobj.lon(Mobj.river_nodes(r))).^2 + (lat - Mobj.lat(Mobj.river_nodes(r))).^2);
    [howclose, idx] = min(dist);
    
    if howclose > dist_thresh
        clim(:, r) = median(climatology, 2);
    else
        % Get the relevant climatology.
        clim(:, r) = climatology(:, idx);
    end
end

% Now we have the climatologies for the relevant river nodes, we need to
% repeat it to fit the length of the Mobj.river_time array.

if ftbverbose
    fprintf('end   : %s \n', subname)
end