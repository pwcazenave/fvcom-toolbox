function Mobj = get_FVCOM_rivers(Mobj, dist_thresh)
% Extract river discharges from the supplied river positions for the FVCOM
% grid in Mobj.
%
% get_FVCOM_rivers(Mobj, rivers, dist_thresh)
%
% DESCRIPTION:
%   For the positioins in fvcom_xy, find the nearest unstructured grid node
%   and extract the river discharge from polcoms_flow. If dist_thresh is
%   specified, the river positions must fall within the specified distance.
%   If multiple rivers are assigned to the same node, their discharges are
%   summed. The resulting river name is generated from the contributing
%   rives, separated by a hyphen.
%
% INPUT:
%   Mobj - MATLAB mesh object containing:
%       * lon, lat - positions for the unstructured grid.
%       * rivers.year - start year of the river data time series.
%       * xy - river positions in lon, lat.
%       * names - list of river names (whose order must match the positions 
%           in xy).
%       * flow - river discharge data (again, order of columns must match
%           the positions in fvcom_xy).
%   dist_thresh - [optional] maximum distance away from a river node beyond
%       which the search for an FVCOM node is abandoned. Units in degrees.
%
% OUTPUT:
%   Mobj.river_flux - volume flux at the nodes within the model domain.
%   Mobj.river_nodes - node IDs for the rivers. At the moment, these are
%       point sources only. Eventually, some rivers may have to be split
%       over several nodes.
%   Mobj.river_names - river names which fall within the model domain. For
%       rivers where the discharge has been summed, the name is compoud,
%       with each contributing name separated by a hyphen (-).
%   Mobj.river_time - time series for the river discharge data
%
% EXAMPLE USAGE:
%   Mobj = get_FVCOM_rivers(Mobj, Mobj.rivers, 0.025)
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2013-03-27 - First version.
%   2013-04-15 - Removed the code to load the positions and discharge data
%   into separate functions so that this function can be purely about
%   finding the closest locations and the summing of discharges (if
%   necessary). The downside is the order of the discharges (columns) must
%   match the position arrays.
%
%==========================================================================

subname = 'get_FVCOM_rivers';

global ftbverbose;
if ftbverbose
    fprintf(['\nbegin : ' subname '\n'])
end

% Check inputs
if ~Mobj.have_lonlat
    error('Require unstructured grid positions in lon/lat format to compare against supplied river positions.')
end

% Separate the inputs into separate arrays.
fvcom_name = Mobj.rivers.names;
fvcom_xy = Mobj.rivers.positions;
polcoms_flow = Mobj.rivers.discharge;

% We have to be careful because POLCOMS has duplicate river names.
%
%   "This has made a lot of people very angry and has been widely regarded
%   as a bad move."
%
% For duplicates, we need, therefore, to work out a way to handle them
% elegantly. We will assume that rivers with the same name are close to one
% another. As such, we'll sum their discharges. 
[~, di] = unique(fvcom_name, 'first');
fv_dupes = 1:length(fvcom_name);
fv_dupes(di) = []; % index of duplicates (does this work with more than two?)

% Iterate through the list of duplicates and clear them out in the names,
% positions and discharge. We have to ensure that the order of the
% deduplicated positions and discharge are the same, otherwise we'll end up
% with the wrong discharge for a given river position, which would be bad.
for i = 1:length(fv_dupes)
    dup_name = fvcom_name(fv_dupes(i));
    didx = strmatch(dup_name, fvcom_name, 'exact');
    
    % Sum the duplicate rivers' data.
    dup_discharge = sum(polcoms_flow(:, didx), 2);
    
    % Remove the original values and put the summed data at the end.
    polcoms_flow(:, didx) = [];
    polcoms_flow = [polcoms_flow, dup_discharge];
    
    % Now remove the duplicates from the FVCOM data.
    fvcom_name{length(fvcom_name) + 1} = fvcom_name{fv_dupes(i)};
    fvcom_name(didx) = [];
    fvcom_xy = [fvcom_xy; fvcom_xy(fv_dupes(i), :)];
    fvcom_xy(didx, :) = [];
end

% Get number of times and rivers from the deduplicated data.
fv_nr = length(fvcom_name);
[pc_nt, ~] = size(polcoms_flow);

clear didx dup_discharge


% Check each location in the FVCOM rivers file against the grid in Mobj and
% for the indices within the dist_thresh, extract the relevant time series
% data.

vc = 0; % valid FVCOM boundary node counter

% We need to find the unstructured grid boundary nodes and exclude the open
% boundary nodes from them. This will be our list of potential candidates
% for the river nodes (i.e. the land coastlin).
[~, ~, ~, bnd] = connectivity([Mobj.lon, Mobj.lat], Mobj.tri);
boundary_nodes = 1:Mobj.nVerts;
boundary_nodes = boundary_nodes(bnd);
coast_nodes = boundary_nodes(~ismember(boundary_nodes, [Mobj.read_obc_nodes{:}]));

fv_obc = nan;
fvcom_names = cell(0);
fv_riv_idx = nan;

for ff = 1:fv_nr
    % Find the open boundary node closest to this river.
    fv_dist = sqrt( ...
        (fvcom_xy(ff, 1) - Mobj.lon(coast_nodes)).^2 + ...
        (fvcom_xy(ff, 2) - Mobj.lat(coast_nodes)).^2);
    [c, idx] = min(fv_dist);
    if c > dist_thresh && dist_thresh ~= -1 % -1 is for no distance check
        if ftbverbose
            fprintf('\tskipping river %s (%f, %f)\n', fvcom_name{ff}, fvcom_xy(ff, 1), fvcom_xy(ff, 2))
        end
        continue
    else
        if ftbverbose
            fprintf('candidate river %s found (%f, %f)... ', fvcom_name{ff}, fvcom_xy(ff, 1), fvcom_xy(ff, 2))
        end
    end

    % Add it to the list of valid rivers
    vc = vc + 1;

    fv_obc(vc) = coast_nodes(idx);

    % We are assuming that the river discharge data array y-dimension is
    % ordered the same as the positions in fvcom_xy. If they are not, then
    % the discharges for the rivers will be incorrect (i.e. you might put
    % the Severn discharge somewhere in the Baltic).
    fvcom_names{vc} = fvcom_name{ff};
    fv_riv_idx(vc) = ff;
    fv_flow(:, vc) = polcoms_flow(:, ff);
    if ftbverbose
        fprintf('added.\n')
    end
end

% Now we've got a list and some of the nodes will be duplicates. Sum the
% discharge values assigned to those nodes.
fv_uniq_obc = unique(fv_obc);

fv_uniq_flow = nan(pc_nt, length(fv_uniq_obc));
fv_uniq_names = cell(length(fv_uniq_obc), 1);

fv_idx = 1:length(fvcom_names);
for nn = 1:length(fv_uniq_obc)
    
    dn = fv_idx(fv_obc == fv_uniq_obc(nn));
    
    fv_uniq_flow(:, nn) = sum(fv_flow(:, dn), 2);
    % Concatenate the river names so we know at least which rivers'
    % discharges have been summed.
    s = fvcom_names(dn);
    s = [sprintf('%s-', s{1:end-1}, s{end})];
    fv_uniq_names{nn} = s(1:end-1); % lose the trailing -.

end

% Assign the relevant arrays to the Mobj.
Mobj.river_nodes = fv_uniq_obc;
Mobj.river_flux = fv_uniq_flow;
Mobj.river_names = fv_uniq_names;

% Create a Modified Julian Day time series starting at January 1st for the
% year in Mobj.rivers.year.
rtimes = datevec( ...
    datenum([Mobj.rivers.year, 1, 1, 0, 0, 0]): ...
    datenum([Mobj.rivers.year, 1, 1, 0, 0, 0]) + pc_nt - 1 ...
    );
Mobj.river_time = nan(pc_nt, 1);
for tt = 1:pc_nt
    Mobj.river_time(tt) = greg2mjulian( ...
        rtimes(tt, 1), rtimes(tt, 2), rtimes(tt, 3), ...
        rtimes(tt, 4), rtimes(tt, 5), rtimes(tt, 6) ...
        );
end

% Figure to check what's going on with identifying river nodes
% figure
% plot(fvcom_xy(:, 1), fvcom_xy(:, 2), 'o', 'MarkerFaceColor', 'b')
% hold on
% plot(Mobj.lon(bnd), Mobj.lat(bnd), 'go', 'MarkerFaceColor', 'g')
% axis('equal', 'tight')
% plot(Mobj.lon(coast_nodes), Mobj.lat(coast_nodes), 'ro')
% plot(Mobj.lon(Mobj.river_nodes), Mobj.lat(Mobj.river_nodes), 'ko', 'MarkerFaceColor', 'k')
% text(Mobj.lon(Mobj.river_nodes) + 0.025, Mobj.lat(Mobj.river_nodes) + 0.025, Mobj.river_names)
% axis([min(Mobj.lon), max(Mobj.lon), min(Mobj.lat), max(Mobj.lat)])
% legend('POLCOMS nodes', 'Grid boundary', 'Land nodes', 'Selected nodes', 'Location', 'NorthOutside', 'Orientation', 'Horizontal')
% legend('BoxOff')


if ftbverbose
    fprintf(['end   : ' subname '\n'])
end
