function Mobj = get_POLCOMS_rivers(Mobj, polcoms_file, polcoms_ij, polcoms_grid, dist_thresh)
% Parse the POLCOMS rivers data file.
%
% get_POLCOMS_rivers(Mobj, polcoms_file, polcoms_ij, polcoms_grid, dist_thresh)
%
% DESCRIPTION:
%   Takes an existing POLCOMS formatted river discharge file and extracts
%   the data which fall within the model domain specified in Mobj.
%
% INPUT:
%   Mobj - MATLAB mesh object containing:
%       * x, y or lon, lat - positions for the unstructured grid.
%   polcoms_file - flow data file(s) from POLCOMS. For multiple files, give
%       in chronological order as a cell array of file names.
%   polcoms_ij - indices in the POLCOMS grid at which each river is
%       located.
%   polcoms_grid - NetCDF file of the POLCOMS grid.
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
%   Mobj = get_POLCOMS_river(Mobj, 'polcoms.flw', 'polcoms.index', ...
%       'polcoms.nc', 0.025)
% 
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2013-03-21 - First version.
%
%==========================================================================

subname = 'get_POLCOMS_rivers';

global ftbverbose;
if ftbverbose
    fprintf(['\nbegin : ' subname '\n'])
end

% Check inputs
if exist(polcoms_grid, 'file') ~= 2
    error('file: %s does not exist or is not readable.', polcoms_grid)
end
if exist(polcoms_ij, 'file') ~= 2
    error('file: %s does not exist or is not readable.', polcoms_ij)
end
if exist(polcoms_file, 'file') ~= 2
    error('file: %s does not exist or is not readable.', polcoms_file)
end
if nargin < 4
    % No distance threshold
    dist_thresh = -1;
end
if ~Mobj.have_lonlat
    error('Require unstructured grid positions in lon/lat format to compare against POLCOMS positions.')
end


% The POLCOMS river file has a pretty straightforward format of a 2D array
% of river along x and time along y. Since it's a simple text file with no
% weird format, we'll just read it in with load.
if iscell(polcoms_file)
    pc_riv = [];
    for rr = 1:length(polcoms_file)
        pc_riv = [pc_riv; load(polcoms_file{rr})];
    end
    clear rr
else
    pc_riv = load(polcoms_file);
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

fidx = fopen(polcoms_ij,'r');
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
            warning('Number of rivers in the index file\n\n\t%s\n\ndoes not match the number of rivers in the data file\n\n\t%s\n', polcoms_ij, polcoms_file)
            warning('Resizing the arrays to match the index file.')
            pc_idx = nan(nridx, 3);
            pc_name = cell(nridx, 1);
            pc_nr = nridx; % reset number of rivers value
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

% Now read in the NetCDF file and grab the real coordinates of those
% positions.
nc = netcdf.open(polcoms_grid, 'NOWRITE');
[~, numvars, ~, ~] = netcdf.inq(nc);

for ii = 1:numvars
    
    [varname, ~, ~, ~] = netcdf.inqVar(nc, ii - 1);
    varid = netcdf.inqVarID(nc, varname);

    if ftbverbose
        fprintf('\tvariable %s... ', varname)
    end
 
    switch varname
        case 'lon'
            pc_lon = netcdf.getVar(nc, varid);
            
        case 'lat'
            pc_lat = netcdf.getVar(nc, varid);

    end
    
    if ftbverbose
        fprintf('done.\n')
    end

end

clear numdims numvars dimnames 

netcdf.close(nc)

% Use the indices from the polcoms_ij file to get the real positions of the
% river nodes.
pc_riv_lonlat = [pc_lon(pc_idx(:, 2), 1), pc_lat(1, pc_idx(:, 3))'];

% Now we have the positions of the rivers, we need to find the most
% appropriate node in the FVCOM grid open boundaries. This will probably
% not work completely because POLCOMS' grid resolution is lower than
% FVCOM's is likely to be, and so the nearest open boundary node might be a
% long way down an estuary, for example. So, we'll do this automatically
% here, but then you'll probably have to move the positions yourself after
% this has finished. Short of some complex searches along the open
% boundaries and some angle-based criteria (sharpest angle = river mouth),
% this is probably the most sensible approach.

vc = 0; % valid FVCOM boundary node counter

% We need to find the coastline nodes and exclude the open boundary nodes
% from them. This will be out list of potential candidates for the river
% nodes.
[~, ~, ~, bnd] = connectivity([Mobj.lon, Mobj.lat], Mobj.tri);
boundary_nodes = 1:Mobj.nVerts;
boundary_nodes = boundary_nodes(bnd);
coast_nodes = boundary_nodes(~ismember(boundary_nodes, [Mobj.read_obc_nodes{:}]));

fv_obc = nan;
fv_names = cell(0);
fv_riv_idx = nan;

for pp = 1:pc_nr
    
    % Find the open boundary node closest to this river
    fv_dist = sqrt( ...
        (pc_riv_lonlat(pp, 1) - Mobj.lon(coast_nodes)).^2 + ...
        (pc_riv_lonlat(pp, 2) - Mobj.lat(coast_nodes)).^2);
    [c, idx] = min(fv_dist);
    if c > dist_thresh && dist_thresh ~= -1 % -1 is for no check
        if ftbverbose
            fprintf('\tskipping river %s (%f, %f)\n', pc_name{pp}, pc_riv_lonlat(pp, 1), pc_riv_lonlat(pp, 2))
        end
        continue
    else
        if ftbverbose
            fprintf('candidate river %s found (%f, %f)... ', pc_name{pp}, pc_riv_lonlat(pp, 1), pc_riv_lonlat(pp, 2))
        end
    end

    % Don't have duplicate river nodes. Is this wise? Probably not because
    % the order in which we stumble upon rivers is arbitrary (alphabetical)
    % and so we might end up with the wrong discharge values. I can't
    % see a more sensible approach than this, though.
    if ~ismember(coast_nodes(idx), fv_obc)
        vc = vc + 1;
        if ftbverbose
            fprintf('added.\n')
        end
        fv_obc(vc) = coast_nodes(idx);
        fv_names{vc} = pc_name{pp};
        fv_riv_idx(vc) = pp;
    else
        if ftbverbose
            fprintf('skipped.\n')
        end
    end
end

% Add to the Mobj.
Mobj.river_nodes = fv_obc;

% Add the relevant time series of river discharge.
Mobj.river_flux = pc_riv(:, fv_riv_idx);

% Time's a bit more interesting because we have approximately daily values
% (the number of rows in pc_riv is approximately number of days in a year,
% plus a few extras). Since the input file name gives us the year, we
% should be able to reconstruct something approaching a date for the
% current time series. As ever, hacking around with input file names is
% bound to break at some point in the future, so warn as such.
warning('Extracting year from the input file name. This might break if your file name doesn''t match what is expected.')
if iscell(polcoms_file)
    [~, fn, ~] = fileparts(polcoms_file{1}); % use the first file name
else
    [~, fn, ~] = fileparts(polcoms_file);
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

Mobj.river_names = fv_names;

% Figure to check what's going on with identifying river nodes
% figure
% plot(pc_riv_lonlat(:,1), pc_riv_lonlat(:,2), 'o', 'MarkerFaceColor', 'b')
% hold on
% plot(Mobj.lon(bnd), Mobj.lat(bnd), 'go', 'MarkerFaceColor', 'g')
% axis('equal', 'tight')
% plot(Mobj.lon(coast_nodes), Mobj.lat(coast_nodes), 'ro')
% plot(Mobj.lon(fv_obc), Mobj.lat(fv_obc), 'ko', 'MarkerFaceColor', 'k')
% axis([min(Mobj.lon), max(Mobj.lon), min(Mobj.lat), max(Mobj.lat)])
% legend('POLCOMS nodes', 'Grid boundary', 'Land nodes', 'Selected nodes', 'Location', 'NorthOutside', 'Orientation', 'Horizontal')
% legend('BoxOff')
