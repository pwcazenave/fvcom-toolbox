function data = get_HYCOM_forcing(Mobj, modelTime, varlist)
% Get mean flow, temperature, salinity, surface elevation and denstiy data
% from HYCOM model outputs through their OPeNDAP server.
%
% data = get_HYCOM_forcing(Mobj, modelTime)
%
% DESCRIPTION:
%   Using OPeNDAP, extract the necessary parameters to create an FVCOM
%   forcing file.
%
% INPUT: 
%   Mobj - MATLAB mesh object with the following fields:
%           - have_lonlat - boolean indicating whether lat/long values are
%           present in Mobj.
%           - lon, lat - longitude and latitudes of the model grid nodes.
%   modelTime - Modified Julian Date start and end times
%   varlist - [optional] cell array of variables to download. Use HYCOM
%       names (e.g. ssh, salinity, temperature, u, v, density). If omitted,
%       downloads salinity, temperature and ssh only.
%
% OUTPUT:
%   data - struct of the data necessary to force FVCOM. These can be
%   interpolated onto the unstructured grid in Mobj using grid2fvcom.m.
%
% The parameters (and the corresponding field names returned) which are
% obtained from the HYCOM data are:
%     - time [MT]
%     - temperature [temperature]
%     - salinity [salinity]
%     - u mean flow component [currently disabled] [u]
%     - v mean flow component [currently disabled] [v]
%     - density [currently disabled] [density]
%     - daily mean sea surface height [ssh]
%
% EXAMPLE USAGE:
%   To download the default parameters (temperature, salinity and ssh):
%
%       modeltime = [55561, 55591]; % time period in Modified Julian Days
%       hycom = get_HYCOM_forcing(Mobj, modeltime);
%
%   To download only sea surface height and density:
%
%       modeltime = [55561, 55591]; % time period in Modified Julian Days
%       hycom = get_HYCOM_forcing(Mobj, modeltime, {'ssh', 'density'})
%
% Author(s)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2013-01-31 First version.
%   2013-08-19 Make some more progress in getting this working. Mainly
%   change the way the dates are handled to form the relevant URL for
%   downloading the data.
%   2013-09-03 More incremetal changes to get the function working. At the
%   moment, I'm ignoring the old OPeNDAP toolbox for reading the data from
%   the HYCOM OPeNDAP server.
%   2013-09-05 It's working! Also add a data.time variable with the time
%   stamps from the HYCOM data.
%   2013-12-02 Add sea surface height to the list of variables that can be
%   downloaded.
%   2013-12-09 Add ability to specify particular variables to download.
%
%==========================================================================

subname = 'get_HYCOM_forcing';

global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

% For checking whether to use the third-party OPeNDAP toolbox or native
% MATLAB tools. OPeNDAP support is included in MATLAB version 7.14 onwards.
% Although this check is made, I haven't actually written the code for the
% third-party toolbox. If you need it, I be pleased if you wrote an
% equivalent to the native version.
v714date = datenum(2012, 3, 1);
currdate = ver('MATLAB');
currdate = datenum(currdate.Date);

if datenum(currdate) < v714date
    error(['This version of MATLAB does not have native OPeNDAP ', ...
        'support and this function relies on that support. If you ', ...
        'require this function, you will have to add the ', ...
        'functionality with the third-party OPeNDAP toolbox.'])
end

% Check if we've been given a cell array of variables and set the varlist
% to that, otherwise default to temperature, salinity and sea surface
% height.
assert(iscell(varlist), 'List of variables to extract must be a cell array: {''var1'', ''var2''}')
if nargin == 2
    varlist = {'temperature', 'salinity', 'ssh'};
end

% Get the extent of the model domain (in spherical).
if ~Mobj.have_lonlat
    error('Need spherical coordinates to extract the forcing data')
else
    % Add a buffer of two grid cells in latitude and two in longitude to
    % make sure the model domain is fully covered by the extracted data.
    [dx, dy] = deal(1/12, 1/12); % HYCOM resolution in degrees
    % West, east, south, north
    extents = [min(Mobj.lon(:)) - (2*dx), max(Mobj.lon(:)) + (2*dx), ...
        min(Mobj.lat(:)) - (2*dy), max(Mobj.lat(:)) + (2*dy)];
end

% Get the URL to use for the first time step.
url = get_url(modelTime(1));

% List of URL suffixes so we can dynamically build the URL for each time
% step.
suffix.MT           = '?MT';
suffix.Longitude    = '?Longitude';
suffix.Latitude     = '?Latitude';
suffix.Depth        = '?Depth';
suffix.temperature  = '?temperature';
suffix.salinity     = '?salinity';
suffix.u            = '?u';
suffix.v            = '?v';
suffix.density      = '?density';
suffix.ssh          = '?ssh';
suffix.X            = '?X';
suffix.Y            = '?Y';

% Set up a struct of the HYCOM data sets in which we're interested.
hycom.MT            = [url, suffix.MT];             % time [1D]
hycom.Longitude     = [url, suffix.Longitude];      % [2D]
hycom.Latitude      = [url, suffix.Latitude];       % [2D]
hycom.Depth         = [url, suffix.Depth];          % water depth [2D]
hycom.temperature   = [url, suffix.temperature];    % [4D]
hycom.salinity      = [url, suffix.salinity];       % [4D]
hycom.ssh           = [url, suffix.ssh];            % sea surface height [3D]
hycom.u             = [url, suffix.u];              % mean flow [4D]
hycom.v             = [url, suffix.v];              % mean flow [4D]
hycom.density       = [url, suffix.density];        % water density [4D]
% hycom.X             = [url, suffix.X];              % crashes MATLAB...
% hycom.Y             = [url, suffix.Y];              % crashes MATLAB...

% Before we go off downloading data, check the variables we've been asked
% for are actually valid HYCOM names.
for vv = 1:length(varlist)
    if iscell(varlist) && ~isfield(hycom, varlist{vv})
        error('Variable %s is not a valid HYCOM variable name.', varlist{vv})
    end
end

data.time = [];

if datenum(currdate) > v714date
    % Use the built in tools to open the remote file.

    % We have to use Longitude here because otherwise MATLAB crashes (for
    % me - version 2012b). Ideally this would be depend on the value of
    % Mobj.nativeCoords.
    ncid = netcdf.open(hycom.Longitude, 'NOWRITE');
    varid = netcdf.inqVarID(ncid, 'Longitude');

    data.X.data = netcdf.getVar(ncid, varid, 'double');

    % Make the longitude values in the range -180 to 180 (to match the
    % model inputs).
    data.X.data = mod(data.X.data, 360);
    data.X.data(data.X.data > 180) = data.X.data(data.X.data > 180) - 360;

    ncid = netcdf.open(hycom.Latitude, 'NOWRITE');
    varid = netcdf.inqVarID(ncid, 'Latitude');

    data.Y.data = netcdf.getVar(ncid, varid, 'double');

    % Create indices of the size of the arrays.
    data.X.idx = repmat(1:size(data.X.data, 1), [size(data.X.data, 2), 1])';
    data.Y.idx = repmat(1:size(data.Y.data, 2), [size(data.Y.data, 1), 1]);
    %data.Y.idx = 1:size(data.Y.data, 2) - 1;

    % Find the indices which cover the model domain then find the extremes
    % (for requesting only a subset from the OPeNDAP server.
    idx = data.X.data > extents(1) & data.X.data < extents(2) & data.Y.data > extents(3) & data.Y.data < extents(4);
    xrange = [min(data.X.idx(idx)), max(data.X.idx(idx))];
    yrange = [min(data.Y.idx(idx)), max(data.Y.idx(idx))];

    data.lon = data.X.data(xrange(1):xrange(2), yrange(1):yrange(2));
    data.lat = data.Y.data(xrange(1):xrange(2), yrange(1):yrange(2));


    % Load the depth data (1D vector).
    if datenum(currdate) > v714date
        ncid = netcdf.open(hycom.Depth);
        varid = netcdf.inqVarID(ncid, 'Depth');

        % HYCOM has fixed depths, so the array which returned here is just
        % a list of those depths. We need them to interpolate the vertical
        % structure onto the FVCOM sigma levels.
        data.Depth.data = netcdf.getVar(ncid, varid, 'double');

    else
        % Third-party toolbox version (unimplemented).
        data_attributes.(fields{aa}) = loaddap('-A', [hycom.(fields{aa})]);
        % Get the data time and convert to Modified Julian Day.
        data.time = loaddap(hycom.time);
    end

    % Get the number of vertical levels.
    nz = length(data.Depth.data);


    % Now find the time indices. HYCOM stores its times as days since
    % 1900-12-31 00:00:00. Naturally this is completely different from
    % everything else ever. We need to iterate through the days we've been
    % given and find the relevant url (with get_url) and then create a cell
    % array of all the times. This assumes the HYCOM data is daily values
    % (which it is).
    times = modelTime(1):modelTime(2);
    nt = length(times);

    % Open the initial connection.
    ncid = netcdf.open(hycom.MT, 'NOWRITE');

    c = 1; % counter for the tmjd cell array.
    for tt = 1:nt
        oldurl = url;
        url = get_url(times(tt));

        % Only reopen the connection if the two URLs differ.
        if ~strcmpi(oldurl, url) || tt == 1

            hycom.MT = [url, suffix.MT];
            ncid = netcdf.open(hycom.MT, 'NOWRITE');
            varid = netcdf.inqVarID(ncid, 'MT');

            % Add the new data to the cell array. This should build up an
            % array of unique time series. We can then use these to query
            % for the time indices for each time step later.
            data.MT.data{c} = netcdf.getVar(ncid, varid, 'double');

            % Convert to Gregorian date and then to Modified Julian Days.
            t{c} = datevec(data.MT.data{c} + datenum(1900, 12, 31, 0, 0, 0));
            tmjd{c} = greg2mjulian(t{c}(:,1), t{c}(:,2), t{c}(:,3), t{c}(:,4), t{c}(:,5), t{c}(:,6));

            c = c + 1;
        end
    end

    % Clear out the full lon/lat arrays.
    data = rmfield(data, {'X', 'Y'});
else
    % I haven't tested this at all.

    data.X.idx = loaddap(hycom.X);
    data.Y.idx = loaddap(hycom.Y);
    xIdx = length(data.X.idx.X) - 1;
    yIdx = length(data.Y.idx.Y) - 1;
    data.lon = loaddap([hycom.lon, sprintf('[%i:1:%i]', 0, 0), sprintf('[%i:1:%i]', 0, xIdx)]);
    data.lat = loaddap([hycom.lat, sprintf('[%i:1:%i]', 0, yIdx), sprintf('[%i:1:%i]', 0, 0)]);
end

netcdf.close(ncid);

% Set the fields to iterate through to be the varlist plus those which are
% essential (time, positions and depth).
% fields = fieldnames(hycom);
fields = ['MT', 'Longitude', 'Latitude', 'Depth', varlist];

for aa = 1:length(fields)

    switch fields{aa}

        case {'Longitude', 'Latitude', 'MT', 'Depth'}
            % Skip these as we've already got the information at the
            % beginning of this function. These are largely time
            % independent data, so no need to get them here as they don't
            % have 4 dimensions, whereas the code after the otherwise
            % assumes 4D data (except for ssh which is 3D).
            continue

        otherwise
            % Store the downloaded data in a struct. Assume the spatial
            % data is identical to that in data.lon and data.lat.
            data.(fields{aa}).data = [];

            if datenum(currdate) > v714date

                ncid = netcdf.open(hycom.(fields{aa}));
                varid = netcdf.inqVarID(ncid, fields{aa});

                % If you don't know what it contains, start by using the
                % 'netcdf.inq' and ncinfo operations:
                %[numdims, numvars, numglobalatts, unlimdimid] = netcdf.inq(ncid);
                %ncid_info = ncinfo(hycom.(fields{aa}));

                % Typically these data are 4D, with dimensions of:
                %   - x (X)
                %   - y (Y)
                %   - depth (Depth)
                %   - time (MT)
                % except in the case of sea surface height, where we lose
                % the depth dimension. For all other variables, we want all
                % depths but only a subset in time and space.

                % Since the HYCOM OPeNDAP server is so spectacularly slow,
                % extract a day's data at a time and stick them together
                % here. If nothing else, this at least means I can give
                % some indication of progress, rather than just wait for
                % something to eventually happen.
                nx = (xrange(2) - xrange(1)) + 1;
                ny = (yrange(2) - yrange(1)) + 1;

                % Preallocate the output so we don't append to an array
                % (which is slow). Sea surface height is 3D only (all the
                % other variables are 4D). So, it needs its own little
                % check all to itself.
                if strcmpi(fields{aa}, 'ssh') == 1
                    was_zeta = true; % set boolean for surface elevation
                    data.(fields{aa}).data = nan(nx, ny, nt);
                else
                    was_zeta = false;
                    data.(fields{aa}).data = nan(nx, ny, nz, nt);
                end

                c = 0; % counter for iterating through tmjd.

                for tt = 1:nt
                    if ftbverbose
                        fprintf('%s: time %i of %i... ', fields{aa}, tt, nt)
                    end

                    % Get the current url value for this time step. This
                    % approach means we can grab data which straddles a
                    % boundary between HYCOM outputs. Only reopen the
                    % connection if the url value has changed. At this
                    % point we also need to get ourselves a new time index
                    % using modelTime and the cell array tmjd.
                    oldurl = url;
                    url = get_url(times(tt));
                    if ~strcmpi(oldurl, url) || tt == 1
                        hycom.(fields{aa}) = [url, suffix.(fields{aa})];

                        ncid = netcdf.open(hycom.(fields{aa}));
                        varid = netcdf.inqVarID(ncid, fields{aa});

                        c = c + 1;
                    end

                    % Find the time index closest to the current model
                    % time.
                    [~, ts] = min(abs(tmjd{c} - times(tt)));

                    if was_zeta
                        % netCDF starts at zero, hence -1.
                        start = [xrange(1), yrange(1), ts] - 1;
                        count = [nx, ny, 1];
                        data.(fields{aa}).data(:, :, tt) = netcdf.getVar(ncid, varid, start, count, 'double');
                    else
                        % netCDF starts at zero, hence -1.
                        start = [xrange(1), yrange(1), 1, ts] - 1;
                        count = [nx, ny, nz, 1];
                        data.(fields{aa}).data(:, :, :, tt) = netcdf.getVar(ncid, varid, start, count, 'double');
                    end

                    % Build an array of the HYCOM times. Only do so once so
                    % we don't end up appending it multiple times.
                    if length(data.time) < nt
                        data.time = [data.time; tmjd{c}(ts)];
                    end

                    if ftbverbose; fprintf('done.\n'); end
                end
            else
                % Third-party toolbox version (unimplemented).

                data_attributes.(fields{aa}) = loaddap('-A', [hycom.(fields{aa})]);
                % Get the data time and convert to Modified Julian Day.
                data.time = loaddap(hycom.time);
            end

            netcdf.close(ncid);
    end
end

if ftbverbose
    fprintf('end   : %s\n', subname)
end

function url = get_url(time)
% Child function to find the relevant URL to use for a given time step.
%
% url = get_url(time);
%
% INPUT:
%   time - Modified Julian Day
%
% OUTPUT:
%   url - string of the approprate URL for the date supplied in time.
%

[t1, t2, t3, t4, t5, t6] = datevec(datestr(now));

if time < greg2mjulian(2008, 09, 16, 0, 0, 0)
    error('Not using the legacy HYCOM model output. Select a start date from September 16th 2008 onwards.')
elseif time >= greg2mjulian(2008, 9, 19, 0, 0, 0) && time < greg2mjulian(2009, 5, 7, 0, 0, 0)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.6';
elseif time >= greg2mjulian(2009, 5, 7, 0, 0, 0) && time < greg2mjulian(2011, 1, 3, 0, 0, 0)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.8';
elseif time >= greg2mjulian(2011, 1, 3, 0, 0, 0) && time < greg2mjulian(2013, 8, 21, 0, 0, 0)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9';
elseif time >= greg2mjulian(2013, 8, 21, 0, 0, 0) && time <= greg2mjulian(t1, t2, t3, t4, t5, t6)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.0';
elseif time > greg2mjulian(t1, t2, t3, t4, t5, t6)
    error('Given date is in the future.')
else
    error('Date is outside of the known spacetime continuum? See help TARDIS.')
end
