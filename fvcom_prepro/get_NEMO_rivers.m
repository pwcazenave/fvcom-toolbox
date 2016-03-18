function Mobj = get_NEMO_rivers(Mobj, dist_thresh, varargin)
% Extract river data from the supplied river positions for the FVCOM
% grid in Mobj from the NEMO data
%
% get_NEMO_rivers(Mobj, dist_thresh)
%
% DESCRIPTION:
%   For the positions in Mobj.rivers.positions, find the nearest
%   unstructured grid node and extract the river discharge from
%   Mobj.rivers.river_flux. The river positions must fall within the
%   specified distance (dist_thresh). If multiple rivers are assigned to
%   the same node, the river with the larger of the set of discharges is
%   used.
%
% INPUT:
%   Mobj - MATLAB mesh object containing:
%       * have_lonlat - boolean to check for spherical coordinates.
%       * lon, lat - positions for the unstructured grid.
%       * tri - triangulation table for the unstructured grid.
%       * nVerts - number of nodes in the grid.
%       * read_obc_nodes - open boundary node IDs.
%       * rivers - river data struct with the following fields:
%           - positions - river positions in lon, lat.
%           - names - list of river names
%           - river_flux - path to the NEMO netCDF data file.
%   dist_thresh - maximum distance away from a river node beyond
%       which the search for an FVCOM node is abandoned. Units in degrees.
%   model_year - [optional] when giving climatology, a year must be
%       specified so that the time series can be anchored in time. The
%       returned time series will be 3 years long centred on the specified
%       year. Discharges will be repeated for the two additional years.
%
% OUTPUT:
%   Mobj.river_flux - volume flux at the nodes within the model domain.
%   Mobj.river_nh4 - ammonia
%   Mobj.river_no3 - notrate
%   Mobj.river_o - oxygen
%   Mobj.river_p - phosphate
%   Mobj.river_sio3 - silicate
%   Mobj.river_dic - dissolved inorganic carbon
%   Mobj.river_bioalk - alkalinity
%   Mobj.river_temp - temperature
%   Mobj.river_salt - salinity
%   Mobj.river_nodes - node IDs for the rivers. At the moment, these are
%       point sources only. Eventually, some rivers may have to be split
%       over several nodes.
%   Mobj.river_names - river names which fall within the model domain. For
%       rivers where the discharge has been summed, the name is compoud,
%       with each contributing name separated by a hyphen (-).
%   Mobj.river_time - time series for the river discharge data
%
% EXAMPLE USAGE:
%   Mobj = get_NEMO_rivers(Mobj, 0.15)
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2016-03-02 - First version based on get_nemo_rivers.m.
%
%==========================================================================

[~, subname] = fileparts(mfilename('fullpath'));

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

% Check inputs
if ~Mobj.have_lonlat
    error(['Require unstructured grid positions in lon/lat format ', ...
        'to compare against supplied river positions.'])
end

yr = [];

% If we have only three arguments, we have to assume we've been given a
% year for the climatology. Otherwise, we need to read the arguments based
% on keyword-value pairs. Really, we want keyword-value pairs all the time,
% so silently work when given three arguments and don't mention it in the
% help. This is going to bite me at some point in the future, I'm sure.
if nargin == 3
    yr = varargin{1};
elseif nargin > 3
    for aa = 1:2:length(varargin)
        switch varargin{aa}
            case 'model_year'
                yr = varargin{aa + 1};
        end
    end
end

if ~isempty(yr) && ~isnumeric(yr)
    error('Trying to do climatology, but don''t have an anchor year. Supply one via the ''model_year'' keyword-value pair.')
end

% Load the NEMO data.
nemo.lon = ncread(Mobj.rivers.river_flux, 'x');
nemo.lat = ncread(Mobj.rivers.river_flux, 'y');
[nemo.LON, nemo.LAT] = meshgrid(nemo.lon, nemo.lat);
nemo.time = ncread(Mobj.rivers.river_flux, 'time_counter');
nemo.flux = ncread(Mobj.rivers.river_flux, 'rorunoff');
nemo.nh4 = ncread(Mobj.rivers.river_flux, 'ronh4');
nemo.no3 = ncread(Mobj.rivers.river_flux, 'rono3');
nemo.o = ncread(Mobj.rivers.river_flux, 'roo');
nemo.p = ncread(Mobj.rivers.river_flux, 'rop');
nemo.sio3 = ncread(Mobj.rivers.river_flux, 'rosio2');
nemo.dic = ncread(Mobj.rivers.river_flux, 'rodic');
nemo.bioalk = ncread(Mobj.rivers.river_flux, 'robioalk');

% Convert units from grams to millimoles where appropriate.
nemo.no3 = (nemo.no3 / 14) * 1000;
nemo.o = (nemo.o / 32) * 1000; % 2 x 16 for O2
nemo.p = (nemo.p / 35.5) * 1000;
nemo.sio3 = (nemo.sio3 / 28) * 1000;

% Flux in NEMO is specified in kg/m^{2}/s. FVCOM wants m^{3}/s. Fix that
% here.
nemo.flux = nemo.flux * 1000;

[~, ~, nt] = size(nemo.flux);

% All the units need to go from per litre to per metre cubed.

names = fieldnames(nemo);
for f = 1:length(names)
    % Skip position data and do nutrient/flux data only.
    switch names{f}
        case {'lon', 'lat', 'LON', 'LAT', 'time'}
            continue
    end
    nemo.(names{f}) = nemo.(names{f}) / 1000;
end

% Now we've got the data, use the flux data to find the indices of the
% rivers in the arrays and extract those as time series in a format
% suitable for writing out with write_FVCOM_river.
mask = sum(nemo.flux, 3) ~= 0;
% Get the indices so we can extract the time series.
[mirow, micol] = ind2sub(size(mask), find(mask == true));
nr = length(mirow);

% Now do all the data.
for n = 1:length(names)
    switch names{n}
        case {'lon', 'lat', 'LON', 'LAT', 'time'}
            continue
    end

    % Preallocate and then fill with the time series for each valid river
    % location.
    nemo.rivers.(names{n}) = nan(nt, nr);
    for r = 1:nr
        nemo.rivers.(names{n})(:, r) = squeeze(nemo.(names{n})(mirow(r), micol(r), :));
    end
end

% Do the positions in the same way (with a loop) to avoid having to figure
% out if using the mask collapses the arrays in teh same way as the find
% did.
nemo.rivers.positions = nan(nr, 2);
nemo.rivers.names = {};
for r = 1:nr
    % For some reason, the output of meshgrid is the wrong way around, so
    % use the row and column indices the wrong way around too.
    nemo.rivers.positions(r, :) = [nemo.LON(micol(r), mirow(r)), ...
        nemo.LAT(micol(r), mirow(r))];
    % We don't have sensible names for the NEMO rivers, so use the position
    % instead. Perhaps there exists a list of the names somewhere which I
    % could use. Yuri has asked Sarah Wakelin to see if she's got that
    % list.
    nemo.rivers.names{r, 1} = sprintf('river_%.6f-%.6f', ...
        nemo.rivers.positions(r, :));
end

% Separate the inputs into separate arrays.
nemo_name = nemo.rivers.names;
nemo_xy = nemo.rivers.positions;

fv_nr = length(nemo_name);

% Check each location in the NEMO positions against the grid in Mobj and
% for the indices within the dist_thresh, load and extract the relevant
% time series data.

vc = 0; % valid FVCOM boundary node counter

% We need to find the unstructured grid boundary nodes and exclude the open
% boundary nodes from them. This will be our list of potential candidates
% for the river nodes (i.e. the land coastline).
[~, ~, ~, bnd] = connectivity([Mobj.lon, Mobj.lat], Mobj.tri);
boundary_nodes = 1:Mobj.nVerts;
boundary_nodes = boundary_nodes(bnd);
coast_nodes = boundary_nodes(~ismember(boundary_nodes, [Mobj.read_obc_nodes{:}]));
tlon = Mobj.lon(coast_nodes);
tlat = Mobj.lat(coast_nodes);

fv_obc = nan;
fv_names = cell(0);

% Initialise the flow array with a 366 day long time series of nans. This
% array will be appended to (unless all rivers are outside the domain).
% Only do this if we're doing climatology (signified by a non-empty year).
skipped = 0;

% Preallocate all the arrays.
for n = 1:length(names)
    switch names{n}
        case {'lon', 'lat', 'LON', 'LAT', 'time'}
            continue
    end
    fv.(names{n}) = nan(nt, nr);
end

for ff = 1:fv_nr
    % Find the coastline node closest to this river. Don't bother with sqrt
    % for the distance threshold, instead just square the distance when
    % doing the comparison. This should increase performance, although
    % probably only marginally.
    fv_dist = (nemo_xy(ff, 1) - tlon).^2 + ...
        (nemo_xy(ff, 2) - tlat).^2;
    [c, idx] = min(fv_dist);
    if c > dist_thresh^2
        skipped = skipped + 1;
        continue
    else
        if ftbverbose
            fprintf('candidate river %s found (%f, %f)... ', nemo_name{ff}, nemo_xy(ff, 1), nemo_xy(ff, 2))
        end
    end

    vc = vc + 1;

    % We need to make sure the element in which this node occurs does not
    % have two land boundaries (otherwise the model just fills up that
    % element because that element will always have a zero velocity).

    % Find the other nodes which are joined to the node we've just found.
    % We don't need the column to get the other nodes in the element, only
    % the row is required.
    [row, ~] = find(Mobj.tri == coast_nodes(idx));

    if length(row) == 1
        % This is a bad node because it is a part of only one element. The
        % rivers need two adjacent elements to work reliably (?). So, we
        % need to repeat the process above until we find a node that's
        % connected to two elements. We'll try the other nodes in the
        % current element before searching the rest of the coastline (which
        % is computationally expensive).

        % Remove the current node index from the list of candidates (i.e.
        % leave only the two other nodes in the element).
        mask = Mobj.tri(row, :) ~= coast_nodes(idx);
        n_tri = Mobj.tri(row, mask);

        % Remove values which aren't coastline values (we don't want to set
        % the river node to an open water node).
        n_tri = intersect(n_tri, coast_nodes);

        % Of the remaining nodes in the element, find the closest one to
        % the original river location (in fvcom_xy).
        [~, n_idx] = sort(sqrt( ...
            (nemo_xy(ff, 1) - Mobj.lon(n_tri)).^2 ...
            + (nemo_xy(ff, 2) - Mobj.lon(n_tri)).^2));

        [row_2, ~] = find(Mobj.tri == n_tri(n_idx(1)));
        if length(n_idx) > 1
            [row_3, ~] = find(Mobj.tri == n_tri(n_idx(2)));
        end
        % Closest first
        if length(row_2) > 1
            idx = find(coast_nodes == n_tri(n_idx(1)));
        % The other one (only if we have more than one node to consider).
        elseif length(n_idx) > 1 && length(row_3) > 1
            idx = find(coast_nodes == n_tri(n_idx(2)));
        % OK, we need to search across all the other coastline nodes.
        else
            % TODO: Implement a search of all the other coastline nodes.
            % My testing indicates that we never get here (at least for the
            % grids I've tested). I'd be interested to see the mesh which
            % does get here...
            continue
        end
        if ftbverbose
            fprintf('alternate node ')
        end
    end

    % Add it to the list of valid rivers
    fv_obc(vc) = coast_nodes(idx);
    fv_names{vc} = nemo.rivers.names{ff};
    
    % Add the current river data to the relevant arrays.
    for n = 1:length(names)
        switch names{n}
            case {'lon', 'lat', 'LON', 'LAT', 'time'}
                continue
        end
        fv.(names{n})(:, vc) = nemo.rivers.(names{n})(:, ff);
    end

    if ftbverbose
        fprintf('added (%f, %f)\n', Mobj.lon(fv_obc(vc)), Mobj.lat(fv_obc(vc)))
    end

end

% Trim the data arrays for the rivers we've extracted since we preallocated
% the array assuming we'd use all the rivers.
for n = 1:length(names)
    switch names{n}
        case {'lon', 'lat', 'LON', 'LAT', 'time'}
            continue
    end
    fv.(names{n}) = fv.(names{n})(:, 1:vc);
end

% Save a list of the field names in the FVCOM river data.
fnames = fieldnames(fv);

% Assign the relevant arrays to the Mobj. Flux is added in the section
% dealing with either climatology or time series data.
Mobj.river_nodes = fv_obc;
Mobj.river_names = fv_names;
Mobj.have_rivers = true;
Mobj.nRivers = length(fv_obc);

% Dump the data into the relevant arrays in Mobj.
for n = 1:length(fnames)
    new = sprintf('river_%s', fnames{n});
    Mobj.(new) = fv.(fnames{n});
end

% Create a Modified Julian Day time series of the NEMO river data. Assume
% it's a climatology, so generate the right time series based on that
% assumption.
checkdate = nemo.time(1:nt) / (60 * 60 * 24);
if isempty(yr) && ~isnumeric(yr)
    error('For climatology, a year must be specified for the time series to be generated.')
end

% Make three years of data starting from the year before the current one.
% Do so accounting for leap years.
daysinyr = [sum(eomday(yr - 1, 1:12)), ...
    sum(eomday(yr, 1:12)), ...
    sum(eomday(yr + 1, 1:12))];
% Offset the checkdate by one to add zero for the first day. Alternative
% would be to specify the day as the end of the previous month (or a day
% value of zero?).
offsetdays = (1:sum(daysinyr)) - 1;
mtime = datevec(datenum(yr - 1, 1, 1, 0, 0, 0) + offsetdays);
Mobj.river_time = greg2mjulian(mtime(:, 1), mtime(:, 2), ...
    mtime(:, 3), mtime(:, 4), mtime(:, 5), mtime(:, 6));

% Repeat the river data for the climatology before adding to the Mobj.
for n = 1:length(fnames)
    new = sprintf('river_%s', fnames{n});
    Mobj.(new) = [fv.(fnames{n})(1:daysinyr(1), :); ...
                  fv.(fnames{n})(1:daysinyr(2), :); ...
                  fv.(fnames{n})(1:daysinyr(3), :)];
end

if ftbverbose
    fprintf('included %d of %d rivers (skipped %d)\n', ...
        fv_nr - skipped, fv_nr, skipped)
    fprintf('end   : %s \n', subname)
end
