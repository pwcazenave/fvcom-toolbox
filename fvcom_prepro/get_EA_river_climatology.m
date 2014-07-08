function Mobj = get_EA_river_climatology(Mobj, ea, dist_thresh)
% Read river temperature climatologies from the Environment Agency river
% temperature data. If no data are found within the threshold specified, a
% mean climatology from the nearest 30 sites is provided instead.
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
%                   - Mobj.river_time - Modified Julian Day array of the
%                   times for the river discharge data (Mobj.river_flux).
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
%   2014-07-08 Think I've fixed the issue with leap years and incorrectly
%   sized output temperature arrays with multiple years.

subname = 'get_EA_river_climatology';

global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

% Load the position (lon/lat), time, climatology and SiteType variables
% only. Not really bothered about the other variables.
nc = netcdf.open(ea, 'NOWRITE');
varid = netcdf.inqVarID(nc, 'climatology');
climatology = netcdf.getVar(nc, varid, 'single');
varid = netcdf.inqVarID(nc, 'time');
time = netcdf.getVar(nc, varid, 'single');
varid = netcdf.inqVarID(nc, 'lon');
lon = netcdf.getVar(nc, varid, 'single');
varid = netcdf.inqVarID(nc, 'lat');
lat = netcdf.getVar(nc, varid, 'single');
varid = netcdf.inqVarID(nc, 'SiteType');
SiteType = netcdf.getVar(nc, varid);
netcdf.close(nc)
clear varid

% Remove any sites which aren't RIVER in SiteType. This is not pretty but
% relatively speedy, so it'll have to do.
good = [];
for i = 1:size(SiteType, 2)
    if strcmp(strtrim(SiteType(:, i)'), 'RIVER')
        good = [good, i];
    end
end

% Clear out the bad sites.
climatology = climatology(:, good);
lon = lon(good);
lat = lat(good);

% Now find the nearest nodes to the river node positions.
nr = length(Mobj.river_nodes);

clim = nan(length(time), nr);

for r = 1:nr
    dist = sqrt((lon - Mobj.lon(Mobj.river_nodes(r))).^2 + (lat - Mobj.lat(Mobj.river_nodes(r))).^2);
    [howclose, idx] = min(dist);

    if howclose > dist_thresh
        % Find the 30 closest sites and use their data instead.
        [~, idx] = sort(dist);
        clim(:, r) = median(climatology(:, idx(1:30)), 2);
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
warning('Don''t know what''s going on with this here. Check the code to find the end day for the river climatology.')

years = unique(yyyy);
ny = length(years);
if ny == 1

    if mod(years, 4) == 0
        endday = (datenum(yyyy(end), mm(end), dd(end), HH(end), MM(end), SS(end)) - ...
            datenum(max(yyyy), 1, 1, 0, 0, 0)) + 1; % add offset of 1 for MATLAB indexing.
    else
        endday = (datenum(yyyy(end), mm(end), dd(end), HH(end), MM(end), SS(end)) - ...
            datenum(max(yyyy), 1, 1, 0, 0, 0));
    end

    % Subset the river climatology for the right days.
    repclim = clim(startday:endday, :);
else
    % Otherwise build up the time series using the number of unique years
    % we have.
    for y = 1:ny
        % Find the number of days in this year and only extract that number
        % from the climatology.
        tidx = 1:length(yyyy);
        tidx(yyyy ~= years(y)) = [];
        if mod(years(y), 4) == 0
            endday = (datenum(yyyy(tidx(end)), mm(tidx(end)), dd(tidx(end)), HH(tidx(end)), MM(tidx(end)), SS(tidx(end))) - ...
                datenum(max(yyyy(tidx)), 1, 1, 0, 0, 0)) + 1; % add offset of 1 for MATLAB indexing.
        else
            endday = (datenum(yyyy(tidx(end)), mm(tidx(end)), dd(tidx(end)), HH(tidx(end)), MM(tidx(end)), SS(tidx(end))) - ...
                datenum(max(yyyy(tidx)), 1, 1, 0, 0, 0));
        end

        nd = sum(eomday(years(y), 1:12));
        if y == 1
            % This is the part year for the first year. Prepend the
            % existing array with the number of days required.
            repclim = clim(startday:end, :);
        elseif y == ny
            repclim = [repclim; clim(1:endday, :)];
        elseif y ~= 1 || y ~= ny
            % We're in the middle years, so just repeat add the clim array
            % to the end of the previous interation's.
            repclim = [repclim; clim];
        end

        % We need to add an extra couple of day's data to the end of the
        % array for this (leap) year.
        if nd == 366
            repclim = [repclim; repclim(end - 1:end, :)];
        end
    end
end

% Add the temperature climatology to Mobj.
Mobj.river_temp = repclim;

if ftbverbose
    fprintf('end   : %s \n', subname)
end