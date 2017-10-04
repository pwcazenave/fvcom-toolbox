function Mobj = get_NEMO_rivers(Mobj, dist_thresh, varargin)
% Extract river data from the supplied river positions for the FVCOM
% grid in Mobj from the NEMO data
%
% get_NEMO_rivers(Mobj, dist_thresh)
%
% DESCRIPTION:
%   For the positions in Mobj.rivers.positions, find the nearest
%   unstructured grid node and extract the river discharge from
%   Mobj.rivers.river_flux. To correct the flux values from the
%   NEMO data, we need the grid size data (in
%   Mobj.rivers.river_coordinates). The river positions must fall
%   within the specified distance (dist_thresh).
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
%           - river_flux - path to the NEMO netCDF data file. Must contain
%           the area of each grid element for conversion of the fluxes.
%   dist_thresh - maximum distance away from a river node beyond
%       which the search for an FVCOM node is abandoned. Units in degrees.
%   The following keyword-argument pairs are also valid:
%   'model_year' - [optional] when giving climatology, a year must be
%       specified so that the time series can be anchored in time. The
%       returned time series will be 3 years long centred on the specified
%       year. Discharges will be repeated for the two additional years.
%   'dump_positions' - [optional] dump the NEMO river positions to the
%       specified text file.
%   'alternate_positions' - [optional] read in a CSV file with alternate
%       positions for the NEMO rivers. Supply a comma separated file with
%       the new values and then the old values as (xnew,ynew,xold,yold).
%       This is useful if you've manually moved the NEMO rivers onto more
%       realistic locations.
%   'remove_baltic' - [optional] remove the Baltic river inputs. Set to
%       true to remove; defaults to false.
%
% OUTPUT:
%   Mobj.river_flux - volume flux at the nodes within the model domain.
%   Mobj.river_nh4 - ammonia
%   Mobj.river_no3 - notrate
%   Mobj.river_o - oxygen
%   Mobj.river_p - phosphate
%   Mobj.river_sio3 - silicate
%   Mobj.river_dic - dissolved inorganic carbon
%   Mobj.river_totalk - total alkalinity
%   Mobj.river_bioalk - bio-alkalinity
%   Mobj.river_nodes - node IDs for the rivers.
%   Mobj.river_names - river names which fall within the model domain. For
%       rivers where the discharge has been summed, the name is compoud,
%       with each contributing name separated by a hyphen (-).
%   Mobj.river_time - Modified Julian Day time series for the river
%       discharge data.
%   Mobj.river_nemo_location - river locations (NEMO positions).
%
% EXAMPLE USAGE:
%   Mobj = get_NEMO_rivers(Mobj, 0.15)
%
%   To extract the NEMO river locations to file:
%
%   Mobj = get_NEMO_rivers(Mobj, 0.15, 'dump_positions', '/tmp/nemo.txt');
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2016-03-02 - First version based on get_nemo_rivers.m.
%   2016-05-03 - Add option to dump NEMO river locations to text file.
%   2016-06-06 - Remove temperature and salinity from the description as
%   they're really the Baltic Sea inputs only. Also read in the grid area
%   from the new format ERSEM file (variable dA). Add total alkalinity
%   to the outputs.
%   2016-08-10 - Add new option to use a file specifying alternative river
%   positions.
%   2017-10-04 - Add new option to remove the Baltic Sea inputs from the
%   river data.
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
dump_positions = false;
alt_positions = false;
drop_baltic = false;
if nargin == 3
    yr = varargin{1};
elseif nargin > 3
    for aa = 1:2:length(varargin)
        switch varargin{aa}
            case 'model_year'
                yr = varargin{aa + 1};
            case 'dump_positions'
                dump_positions = true;
                position_file = varargin{aa + 1};
            case 'alternate_positions'
                alt_positions = true;
                alternate_file = varargin{aa + 1};
            case 'remove_baltic'
                drop_baltic = true;
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
nemo.alt = ncread(Mobj.rivers.river_flux, 'rototalk');
nemo.bioalk = ncread(Mobj.rivers.river_flux, 'robioalk');
% nemo.temp = ncread(Mobj.rivers.river_flux, 'rotemper');
% nemo.salt = ncread(Mobj.rivers.river_flux, 'rosaline');

if drop_baltic
    % Set the data at the Baltic indices to zero before we get too carried
    % away. Find the indices based on positions rather than hard-coding
    % them so we can still do this for newer NEMO river inputs (assuming
    % they still use the same approach for the Baltic).
    baltic.lon = [10.7777, 12.5555];
    baltic.lat = [55.5998, 56.1331];
    names = fieldnames(nemo);
    for pp = 1:length(baltic.lon)
        [~, lon_idx] = min(abs(nemo.lon - baltic.lon(pp)));
        [~, lat_idx] = min(abs(nemo.lat - baltic.lat(pp)));
        for n = 1:length(names)
            switch names{n}
                case {'lon', 'lat', 'LON', 'LAT', 'time'}
                    continue
            end

            % Drop the Baltic data (replace with zeros to match the other
            % non-river data in the netCDF).
            nemo.(names{n})(lon_idx, lat_idx, :) = 0;
        end
    end
    clearvars lon_idx lat_idx baltic names n
end

% Now get the NEMO grid data.
% nemo.e1t = ncread(Mobj.rivers.river_coordinates, 'e1t');
% nemo.e2t = ncread(Mobj.rivers.river_coordinates, 'e2t');
% Calculate the area for all elements in the grid.
% nemo.area = nemo.e1t .* nemo.e2t;
% Just use the new area variable instead of needing two files.
nemo.area = ncread(Mobj.rivers.river_flux, 'dA');

% NEMO does conversions to ERSEM units internally. Whilst this is easy in
% some ways, it's not particularly transparent. So, instead, we'll do all
% the conversions up front and then the data that get loaded into ERSEM are
% already in the correct units. To summarise those conversions, we have:
%
%   nutrients (nh4, no3, o, p, sio3) from grams/l to mmol/m^3
%   flux from kg/m^{2}/s to m^{3}/s (divide by freshwater density)
%   dic - no change whatsoever.

[~, ~, nt] = size(nemo.flux);

% Flux in NEMO is specified in kg/m^{2}/s. FVCOM wants m^{3}/s. Divide by
% freshwater density to get m/s and then multiply by the area of each
% element to get flux.
nemo.flux = nemo.flux / 1000;
% Now multiply by the relevant area to (finally!) get to m^{3}/s.
nemo.flux = nemo.flux .* repmat(nemo.area, 1, 1, nt);
% Set zero values to a very small number instead.
tmp = nemo.flux;
tmp(tmp==0) = 1E-8;

% Convert units from grams to millimoles where appropriate.
nemo.nh4 = (nemo.nh4 / 14) *1000 ./  tmp; %g/s to mmol/m3
nemo.no3 = (nemo.no3 / 14) *1000 ./  tmp;%g/s to mmol/m3
nemo.o = (nemo.o / 16) *1000 ./ tmp; % Nemo oxygen concentrations are for O rather than O2
nemo.p = (nemo.p / 35.5)*1000 ./ tmp;%g/s to mmol/m3
nemo.sio3 = (nemo.sio3 / 28) *1000./ tmp;%g/2 to mmol/m3
nemo.bioalk = nemo.bioalk./ tmp / 1000; % bioalk is in umol/s need umol/kg
nemo.dic = nemo.dic./12./ tmp *1000; % dic is in gC/s need mmol/m3
% total alkalinity is already in umol/Kg as expected by ERSEM.
clear tmp

% Now we've got the data, use the flux data to find the indices of the
% rivers in the arrays and extract those as time series in a format
% suitable for writing out with write_FVCOM_river.
mask = sum(nemo.flux, 3) ~= 0;
% Get the indices so we can extract the time series.
[mirow, micol] = ind2sub(size(mask), find(mask == true));
nr = length(mirow);

% Now do all the data.
names = fieldnames(nemo);
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

if alt_positions
    % We've been given a file with alternate river positions in. Use that
    % for the river locations instead.
    f = fopen(alternate_file, 'r');
    assert(f > 1, 'Error opening river locations file (%s)', alternate_file)
    new_positions = textscan(f, '%f%f%f%f%[^\n\r]', ...
        'Delimiter', ',', ...
        'HeaderLines', 1, ...
        'ReturnOnError', false);
    new_x = new_positions{1};
    new_y = new_positions{2};
    original_x = new_positions{3};
    original_y = new_positions{4};
    % Although in principle the "old" positions should be identical, odd
    % little precision issues can creep in and cause all sorts of problems.
    % So, we'll instead search for the nearest location and use that
    % instead.
    for ri = 1:length(original_x)
        xdiffs = original_x(ri) - nemo.rivers.positions(:, 1);
        ydiffs = original_y(ri) - nemo.rivers.positions(:, 2);
        [~, index] = min(sqrt(xdiffs.^2 - ydiffs.^2));
        % Now we know which river we're updating, just replace the NEMO
        % positions read in from netCDF with our CSV versions.
        nemo.rivers.positions(index, :) = [new_x(ri), new_y(ri)];
    end
    clear xdiffs ydiffs ri index
end

% Separate the inputs into separate arrays.
nemo_name = nemo.rivers.names;
nemo_xy = nemo.rivers.positions;

if dump_positions
    % Add a header for GIS
    dlmwrite(position_file, 'lonDD,latDD', 'delimiter', '');
    dlmwrite(position_file, nemo_xy, 'precision', '%0.6f', '-append');
end

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
fv_location = [nan, nan];
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
            + (nemo_xy(ff, 2) - Mobj.lat(n_tri)).^2));

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
    fv_location(vc, :) = nemo.rivers.positions(ff, :);

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
Mobj.river_names = fv_names';
Mobj.river_nemo_location = fv_location;
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

    % NEMO is insane. Old versions of the rivers data had enough days for
    % leap and non-leap years; new versions don't. So, let's pad the data
    % by a day in either case and that should work for data with both 365
    % and 366 times. In the case where we're using data with 366 values
    % already, this shouldn't make any difference as we're explicitly
    % indexing to the days in the year using daysinyr anyway.
    fv.(fnames{n}) = [fv.(fnames{n}); fv.(fnames{n})(end, :)];

    Mobj.(new) = [fv.(fnames{n})(1:daysinyr(1), :); ...
                  fv.(fnames{n})(1:daysinyr(2), :); ...
                  fv.(fnames{n})(1:daysinyr(3), :)];
end
Mobj.rivers_orig = nemo.rivers;
if ftbverbose
    fprintf('included %d of %d rivers (skipped %d)\n', ...
        fv_nr - skipped, fv_nr, skipped)
    fprintf('end   : %s \n', subname)
end
