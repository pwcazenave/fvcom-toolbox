function Mobj = get_FVCOM_rivers(Mobj, fvcom_xy, polcoms_flow, polcoms_ij, dist_thresh)
% Extract positions from the FVCOM-adapted POLCOMS rivers files.
%
% get_FVCOM_rivers(Mobj, fvcom_xy, polcoms_flow, polcoms_ij, dist_thresh)
%
% DESCRIPTION:
%   Takes a list of positions and names in ASCII format (fvcom_xy) and
%   finds the nearest node to that position in the grid in Mobj. Discharge
%   values are extracted based on the station name from the POLCOMS
%   formatted flow files. If a two (or more) rivers are assigned to the
%   same node in the FVCOM grid, the discharge values are summed at that
%   node.
%
% INPUT:
%   Mobj - MATLAB mesh object containing:
%       * lon, lat - positions for the unstructured grid.
%   fvcom_xy - ASCII list of lon, lat, name (comma separated).
%   polcoms_flow - flow data file(s) from POLCOMS. For multiple files, give
%       in chronological order as a cell array of file names.
%   polcoms_ij - indices in the POLCOMS grid at which each river is
%       located. The order of these locations must match the order in the
%       flow file as it will be used to determine which time series in
%       polcoms_flow belongs with which station name.
%   dist_thresh - [optional] Maximum distance away from a POLCOMS river
%       node beyond which the search for an FVCOM node is abandoned. Units
%       in degrees.
% 
% OUTPUT:
%   Mobj.river_flux - volume flux at the nodes within the model domain.
%   Mobj.river_nodes - node IDs for the rivers. At the moment, these are
%       point sources only. Eventually, some rivers may have to be split
%       over several nodes.
%   Mobj.river_time - time series for the river discharge data
%
% EXAMPLE USAGE:
%   Mobj = get_FVCOM_rivers(Mobj, 'fvcom_xy.csv', 'polcoms.flw', 0.025)
% 
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2013-03-27 - First version.
%
%==========================================================================

subname = 'get_FVCOM_rivers';

global ftbverbose;
if ftbverbose
    fprintf(['\nbegin : ' subname '\n'])
end

% Check inputs
if exist(fvcom_xy, 'file') ~= 2
    error('file: %s does not exist or is not readable.', fvcom_xy)
end
if exist(polcoms_flow, 'file') ~= 2
    error('file: %s does not exist or is not readable.', polcoms_flow)
end
if nargin < 5
    % No distance threshold
    dist_thresh = -1;
end
if ~Mobj.have_lonlat
    error('Require unstructured grid positions in lon/lat format to compare against POLCOMS positions.')
end

% The POLCOMS river file has a pretty straightforward format of a 2D array
% of river along x and time along y. Since it's a simple text file with no
% weird format, we'll just read it in with load.
if iscell(polcoms_flow)
    pc_riv = [];
    for rr = 1:length(polcoms_flow)
        pc_riv = [pc_riv; load(polcoms_flow{rr})];
    end
    clear rr
else
    pc_riv = load(polcoms_flow);
end
[pc_nt, pc_nr] = size(pc_riv);

% Get the positions for each river from the NetCDF grid and the index file.
% Format is:
% 
% n
%   id1 i j Name1
%   id2 i j Name2
%   ...
%   idn i j Namen

fidx = fopen(polcoms_ij, 'r');
if fidx < 0
    error('file: %s does not exist', polcoms_ij);
end

% Output arrays
pc_idx = nan(pc_nr, 3);
pc_name = cell(pc_nr, 1);

c = 0; % line counter
while ~feof(fidx)
    line = fgetl(fidx);
    if isempty(line) || ~ischar(line)
        continue
    else
        c = c + 1;
    end
    
    if c == 1
        % First (valid) line should be number of rivers
        nridx = str2double(strtrim(line));
        if nridx ~= pc_nr
            warning('Number of rivers in the index file\n\n\t%s\n\ndoes not match the number of rivers in the data file\n\n\t%s\n\n', polcoms_ij, polcoms_flow)
            nridx = pc_nr; % reset number of rivers value
            wrongSize = true;
        else
            wrongSize = false;
        end
    else
        % We're in the data
        S = regexpi(strtrim(line), ' +', 'split');
        pc_idx(c - 1, :) = [str2double(S{1}), str2double(S{2}), str2double(S{3})];
        pc_name{c - 1} = S{end};
    end
end
clear S line c

fclose(fidx);

if wrongSize
    warning('Truncating the POLCOMS name and index arrays to match the data file.')
    pc_name = pc_name(1:pc_nr);
    pc_idx = pc_idx(1:pc_nr, :);
end
% Now get the locations from the fvcom_xy CSV file.
fid = fopen(fvcom_xy, 'r');
if fid < 0
    error('Unable to open positions file %s.', fvcom_xy)
end
fv_xy = textscan(fid, '%f %f %s', 'delimiter', ',', 'HeaderLines', 1);

fclose(fid);

fv_pos = [fv_xy{1}, fv_xy{2}];
fv_name = fv_xy{3};

fv_nr = length(fv_name);

clear fv_xy

if fv_nr ~= pc_nr
    error('The list of FVCOM river positions does not match the number of POLCOMS river time series.')
end


% Find the relevant index for this FVCOM river in the POLCOMS river name
% index (pc_name). We have to be careful here because POLCOMS has duplicate
% river names.
%
%   "This has made a lot of people very angry and has been widely regarded
%   as a bad move."
%
% For duplicates, we need therefore to work out a way to handle them
% elegantly. We will assume that rivers with the same name are close to one
% another and therefore we'll sum their discharges. That way, we can
% dispense with complicated code which looks at whether we've got a
% duplicate and tries to find the existing discharges etc. Instead, we just
% assume that all matches are unique.
[~, di] = unique(fv_name, 'first');
fv_dupes = 1:length(fv_name);
fv_dupes(di) = []; % index of duplicates (does this work with more than two?)

% Iterate through the list of duplicates and pull out the index of the
% other values. Then, sum the two new data sets and delete the original
% ones.
for i = 1:length(fv_dupes)
    dup_name = fv_name(fv_dupes(i));
    didx = strmatch(dup_name, fv_name, 'exact');
    
    % Sum the duplicate rivers' data.
    dup_discharge = sum(pc_riv(:, didx), 2);
    
    % Remove the original values and put the summed data at the end.
    pc_riv(:, didx) = [];
    pc_riv = [pc_riv, dup_discharge];
    
    % Update the list of names to eliminate the duplicates
    pc_name{length(pc_name) + 1} = pc_name{fv_dupes(i)};
    pc_name(didx) = [];
    
    % Clear out the indices too (not sure we're using them here, but best
    % to be thorough). Use the first index even if the positions aren't the
    % same.
    pc_idx = [pc_idx; pc_idx(fv_dupes(1), :)];
    pc_idx(didx, :) = [];
    
    % Now remove the duplicates from the FVCOM data.
    fv_name{length(fv_name) + 1} = fv_name{fv_dupes(1)};
    fv_name(didx) = [];
    fv_pos = [fv_pos; fv_pos(fv_dupes(1), :)];
    fv_pos(didx, :) = [];

end

% Reset the counts to match the deduplicated data.
fv_nr = length(fv_name);
[pc_nt, pc_nr] = size(pc_riv);

clear dup_name didx dup_discharge


% Check each location in the FVCOM rivers file against the grid in Mobj and
% for the indices within the dist_thresh, extract the relevant time series
% data.

vc = 0; % valid FVCOM boundary node counter

% We need to find the coastline nodes and exclude the open boundary nodes
% from them. This will be our list of potential candidates for the river
% nodes.
[~, ~, ~, bnd] = connectivity([Mobj.lon, Mobj.lat], Mobj.tri);
boundary_nodes = 1:Mobj.nVerts;
boundary_nodes = boundary_nodes(bnd);
coast_nodes = boundary_nodes(~ismember(boundary_nodes, [Mobj.read_obc_nodes{:}]));

fv_obc = nan;
fv_names = cell(0);
fv_riv_idx = nan;

% Create an index to find the POLCOMS column for each of the FVCOM rivers.
fv_pc = 1:fv_nr;

for ff = 1:fv_nr
    % Find the open boundary node closest to this river
    fv_dist = sqrt( ...
        (fv_pos(ff, 1) - Mobj.lon(coast_nodes)).^2 + ...
        (fv_pos(ff, 2) - Mobj.lat(coast_nodes)).^2);
    [c, idx] = min(fv_dist);
    if c > dist_thresh && dist_thresh ~= -1 % -1 is for no check
        if ftbverbose
            fprintf('\tskipping river %s (%f, %f)\n', fv_name{ff}, fv_pos(ff, 1), fv_pos(ff, 2))
        end
        continue
    else
        if ftbverbose
            fprintf('candidate river %s found (%f, %f)... ', fv_name{ff}, fv_pos(ff, 1), fv_pos(ff, 2))
        end
    end

    % Add it to the list of valid rivers
    vc = vc + 1;

    fv_obc(vc) = coast_nodes(idx);

    % We need to find the index of the name in the POLCOMS data to make
    % sure we extract the right discharge data (this approach means
    % that the FVCOM data and the POLCOMS data can be in completely
    % different orders and it won't matter (he says...)).
    pc_tmp_idx = strmatch(fv_name{ff}, pc_name, 'exact');
    if length(pc_tmp_idx) > 1
        % Shouldn't get here because we've deduplicated the data
        % sets...
        error('Multiple rivers of name %s.', fv_name{ff})
    end
    fv_names{vc} = pc_name{pc_tmp_idx};
    fv_riv_idx(vc) = pc_tmp_idx;
    fv_flow(:, vc) = pc_riv(:, pc_tmp_idx);
    if ftbverbose
        fprintf('added.\n')
    end
end

% Now we've got a list and some of the nodes will be duplicates. Sum the
% discharge values assigned to those nodes.

fv_uniq_obc = unique(fv_obc);

fv_uniq_flow = nan(pc_nt, length(fv_uniq_obc));
fv_uniq_names = cell(length(fv_uniq_obc), 1);

fv_idx = 1:length(fv_names);
for nn = 1:length(fv_uniq_obc)
    
    dn = fv_idx(fv_obc == fv_uniq_obc(nn));
    
    fv_uniq_flow(:, nn) = sum(fv_flow(:, dn), 2);
    % Concatenate the river names so we know at least which rivers'
    % discharges have been summed.
    s = fv_names(dn);
    s = [sprintf('%s-', s{1:end-1}, s{end})];
    fv_uniq_names{nn} = s(1:end-1); % lose the trailing -.

end

% Assign the relevant arrays to the Mobj.
Mobj.river_nodes = fv_uniq_obc;
Mobj.river_flux = fv_uniq_flow;
Mobj.river_names = fv_uniq_names;

% Time's a bit more interesting because we have approximately daily values
% (the number of rows in pc_riv is approximately number of days in a year,
% plus a few extras). Since the input file name gives us the year, we
% should be able to reconstruct something approaching a date for the
% current time series. As ever, hacking around with input file names is
% bound to break at some point in the future, so warn as such.
warning('Extracting year from the input file name. This might break if your file name doesn''t match what is expected.')
if iscell(polcoms_flow)
    [~, fn, ~] = fileparts(polcoms_flow{1}); % use the first file name
else
    [~, fn, ~] = fileparts(polcoms_flow);
end
% Assume last four characters are the year.
ryear = str2double(fn(end-3:end));

% Create a Modified Julian Day time series starting at January 1st, ryear.
rtimes = datevec(datenum([ryear, 1, 1, 0, 0, 0]):datenum([ryear, 1, 1, 0, 0, 0]) + pc_nt - 1);
Mobj.river_time = nan(pc_nt, 1);
for tt = 1:pc_nt
    Mobj.river_time(tt) = greg2mjulian( ...
        rtimes(tt, 1), rtimes(tt, 2), rtimes(tt, 3), ...
        rtimes(tt, 4), rtimes(tt, 5), rtimes(tt, 6) ...
        );
end

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end
