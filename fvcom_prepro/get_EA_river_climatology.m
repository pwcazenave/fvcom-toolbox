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
% repeat it to fit the length of the Mobj.river_time array. Since the
% climatology data are for a year only, we need to find the correct index
% for the start and end of the Mobj.river_time array so that we don't put
% January temperature in July, for example.
[yyyy, mm, dd, HH, MM, SS] = mjulian2greg(Mobj.river_time);
startday = (datenum(yyyy(1), mm(1), dd(1), HH(1), MM(1), SS(1)) - ...
    datenum(min(yyyy), 1, 1, 0, 0, 0)) + 1; % add offset of 1 for MATLAB indexing.
endday = (datenum(yyyy(end), mm(end), dd(end), HH(end), MM(end), SS(end)) - ...
    datenum(max(yyyy), 1, 1, 0, 0, 0));

years = unique(yyyy);
ny = length(years);
if ny == 1
    % Subset the river climatology for the right days.
    repclim = clim(startday:endday, :);
else
    % Otherwise build up the time series using the number of unique years
    % we have.
    for y = 1:ny
        % Find the number of days in this year and only extract that number
        % from the climatology.
        nd = sum(eomday(years(y), 1:12));
        if y == 1
            % This is the part year for the first year. Prepend the existing
            % array with the number of days required.
            repclim = clim(startday:end, :);
        elseif y == ny
            repclim = [repclim; clim(1:endday, :)];
        elseif y ~= 1 || y ~= ny
            % We're in the middle years, so just repeat add the clim array to
            % the end of the previous interation's.
            repclim = [repclim; clim];
        end

        % We need to add an extra day's data to the end of the array for this
        % year.
        if nd == 366
            repclim = [repclim; repclim(end, :)];
        end
    end
end

% Add the temperature climatology to Mobj.
Mobj.river_temp = repclim;

if ftbverbose
    fprintf('end   : %s \n', subname)
end