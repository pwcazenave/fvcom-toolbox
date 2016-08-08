function data = get_HYCOM_forcing(Mobj, modelTime, varargin)
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
%       names (e.g. ssh, salinity, temperature, u, v). If omitted,
%       downloads salinity, temperature and ssh only.
%   three - [optional] boolean. Set to true to download the 3 hourly data
%       for the period 1992/10/2 to 2008/09/18 (inclusive). Defaults to
%       false and downloads only the daily data.
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
%     - u mean flow component [u]
%     - v mean flow component [v]
%     - daily mean sea surface height [ssh]
%
% EXAMPLE USAGE:
%   To download the default parameters (temperature, salinity and ssh):
%
%       modeltime = [55561, 55591]; % time period in Modified Julian Days
%       hycom = get_HYCOM_forcing(Mobj, modeltime);
%
%   To download only sea surface height:
%
%       modeltime = [55561, 55591]; % time period in Modified Julian Days
%       hycom = get_HYCOM_forcing(Mobj, modeltime, {'ssh'})
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
%   2013-12-12 Fix the handling of the variable input list of field names.
%   2015-05-21 Add support for the Global Reanalysis data which extends
%   coverage back to 1992 (previously limited to 2008 with the Global
%   Analysis data). The Global Analysis data is used from 2008 onwards even
%   though the reanalysis exists up to 2012.
%   2016-01-04 Add support for the three hourly output data for the 19.0
%   and 19.1 experiments.
%   2016-07-06 Add new data sets from the HYCOM server (post-2014).
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
% height. For the three option, assume we want daily data unless three is
% true, in which case set the threehourly parameter to be the string
% required for the THREDDS URL (/3hrly), otherwise default to daily data.
threehourly = '';
if nargin == 2 || nargin > 3
    varlist = {'temperature', 'salinity', 'ssh'};
end
if nargin > 2
    % To maintain backwards compatibility, assume if we've been given a
    % cell array as the third option, that's the variables to download.
    % Otherwise, iterate through the argument as key-parameter pairs.
    if nargin == 3
        if iscell(varargin{1})
            varlist = varargin{1};
            assert(iscell(varargin{1}), [...
                'List of variables to extract must', ...
                ' be a cell array: {''var1'', ''var2''}'])
        end
    else
        for v = 1:2:length(varargin)
            switch varargin{v}
                case 'three'
                    if varargin{v + 1}
                        threehourly = '/3hrly';
                    end
                    if modelTime(1) >= greg2mjulian(2008, 09, 19, 0, 0, 0)
                        error(['The three hourly output is not ', ...
                                'configured for the period beyond ', ...
                                '2008/09/18 in this function. Please ', ...
                                'disable three hourly downloads, or ', ...
                                'edit this function to suit your needs.'])
                    end
            end
        end
    end
end

% Get the extent of the model domain (in spherical).
if ~Mobj.have_lonlat
    error('Need spherical coordinates to extract the forcing data')
else
    % Add a buffer of two grid cells in latitude and two in longitude to
    % make sure the model domain is fully covered by the extracted data.
    [dx, dy] = deal(1/12, 1/12); % HYCOM resolution in degrees
    % West, east, south, north
    extents = [min(Mobj.lon(:)) - (2 * dx), max(Mobj.lon(:)) + (2 * dx), ...
        min(Mobj.lat(:)) - (2 * dy), max(Mobj.lat(:)) + (2 * dy)];
end

% List of URL suffixes so we can dynamically build the URL for each time
% step.
suffix.MT           = {'time', 'MT'};
suffix.Longitude    = {'lon', 'Longitude'};
suffix.Latitude     = {'lat', 'Latitude'};
suffix.Depth        = {'depth', 'Depth'};
suffix.temperature  = {'water_temp', 'temperature'};
suffix.salinity     = {'salinity', 'salinity'};
suffix.u            = {'water_u', 'u'};
suffix.v            = {'water_v', 'v'};
suffix.ssh          = {'surf_el', 'ssh'};

% Get the URL to use for the first time step.
url = get_url(modelTime(1), threehourly);

if modelTime(1) < greg2mjulian(2008, 09, 19, 0, 0, 0)
    name_index = 1;
elseif modelTime(1) >= greg2mjulian(2008, 09, 19, 0, 0, 0)
    name_index = 2;
end
hycom.MT            = [url, suffix.MT{name_index}];          % time [1D]
hycom.Longitude     = [url, suffix.Longitude{name_index}];   % [2D]
hycom.Latitude      = [url, suffix.Latitude{name_index}];    % [2D]
hycom.Depth         = [url, suffix.Depth{name_index}];       % water depth [2D]
hycom.temperature   = [url, suffix.temperature{name_index}]; % [4D]
hycom.salinity      = [url, suffix.salinity{name_index}];    % [4D]
hycom.ssh           = [url, suffix.ssh{name_index}];         % sea surface height [3D]
hycom.u             = [url, suffix.u{name_index}];           % mean flow [4D]
hycom.v             = [url, suffix.v{name_index}];           % mean flow [4D]

% Load the depth data (1D vector).
ncid = netcdf.open(hycom.Depth, 'NOWRITE');
try
    varid = netcdf.inqVarID(ncid, 'Depth');
catch
    varid = netcdf.inqVarID(ncid, 'depth');
end

% HYCOM has fixed depths, so the array which returned here is just
% a list of those depths. We need them to interpolate the vertical
% structure onto the FVCOM sigma levels.
data.Depth.data = netcdf.getVar(ncid, varid, 'double');

netcdf.close(ncid)

% Get the number of vertical levels.
nz = length(data.Depth.data);

% Set a time increment is we've got the three-hourly flag passed to the
% function. At some point, this may have to be a string comparison such
% that we can use different outputs (if the HYCOM folks add them), but for
% now, we just assume non-empty strings means three-hourly data. This will
% also break horribly if we request data from both daily and three-hourly
% data sets.
if ~isempty(threehourly)
    time_increment = 3/24;
else
    time_increment = 1;
end
times = modelTime(1):time_increment:modelTime(2);
nt = length(times);

% Before we go off downloading data, check the variables we've been asked
% for are actually valid HYCOM names.
for vv = 1:length(varlist)
    if iscell(varlist) && ~isfield(hycom, varlist{vv})
        error('Variable %s is not a valid HYCOM variable name.', varlist{vv})
    end
end

% Initialise the value of url so we can compare it each loop and only
% reopen the connection to the server if we've crossed into a different
% data set.
url = get_url(modelTime(1), threehourly);

data.time = [];
c = 1; % counter for the tmjd cell array.
for tt = 1:nt

    % So we can use either the Reanalysis (pre-2008) or Analysis
    % (post-2008) data, we need to build the request struct based on the
    % current time.

    % Set up a struct of the HYCOM data sets in which we're interested.
    if times(tt) < greg2mjulian(2008, 09, 19, 0, 0, 0)
        name_index = 1;
    elseif times(tt) >= greg2mjulian(2008, 09, 19, 0, 0, 0)
        name_index = 2;
    end
    hycom.MT            = [url, suffix.MT{name_index}];          % time [1D]
    hycom.Longitude     = [url, suffix.Longitude{name_index}];   % [2D]
    hycom.Latitude      = [url, suffix.Latitude{name_index}];    % [2D]
    hycom.Depth         = [url, suffix.Depth{name_index}];       % water depth [2D]
    hycom.temperature   = [url, suffix.temperature{name_index}]; % [4D]
    hycom.salinity      = [url, suffix.salinity{name_index}];    % [4D]
    hycom.ssh           = [url, suffix.ssh{name_index}];         % sea surface height [3D]
    hycom.u             = [url, suffix.u{name_index}];           % mean flow [4D]
    hycom.v             = [url, suffix.v{name_index}];           % mean flow [4D]

    oldurl = url;
    url = get_url(times(tt), threehourly);
    % Only reopen the connection if the two URLs differ.
    if ~strcmpi(oldurl, url) || tt == 1
        if times(tt) < greg2mjulian(2008, 09, 19, 0, 0, 0)
            hycom.MT = [url, suffix.MT{1}];
        elseif times(tt) >= greg2mjulian(2008, 09, 19, 0, 0, 0)
            hycom.MT = [url, suffix.MT{2}];
        end
        ncid = netcdf.open(hycom.MT, 'NOWRITE');
        try
            varid = netcdf.inqVarID(ncid, 'MT');
        catch
            varid = netcdf.inqVarID(ncid, 'time');
        end

        % Add the new data to the cell array. This should build up an
        % array of unique time series. We can then use these to query
        % for the time indices for each time step later.
        data.MT.data{c} = netcdf.getVar(ncid, varid, 'double');

        netcdf.close(ncid)

        % Convert to Gregorian date and then to Modified Julian Days. The
        % Global Reanalysis stores time as hours since 2000-01-01, the
        % Global Analysis as days since 1900-12-31.
        if times(tt) < greg2mjulian(2008, 09, 19, 0, 0, 0)
            t{c} = datevec((data.MT.data{c} / 24) + datenum(2000, 1, 1, 0, 0, 0));
        elseif times(tt) >= greg2mjulian(2008, 09, 19, 0, 0, 0)
            t{c} = datevec(data.MT.data{c} + datenum(1900, 12, 31, 0, 0, 0));
        end
        tmjd{c} = greg2mjulian(t{c}(:,1), t{c}(:,2), t{c}(:,3), t{c}(:,4), t{c}(:,5), t{c}(:,6));

        c = c + 1;
    end
end

% Open the relevant spatial variables on the remote server and download the
% spatial data required to subset the HYCOM data.
ncid = netcdf.open(hycom.Longitude, 'NOWRITE');
try
    varid = netcdf.inqVarID(ncid, suffix.Longitude{1});
catch
    varid = netcdf.inqVarID(ncid, suffix.Longitude{2});
end

data.X.data = netcdf.getVar(ncid, varid, 'double');

netcdf.close(ncid)

% Make the longitude values in the range -180 to 180 (to match the
% model inputs).
data.X.data = mod(data.X.data, 360);
data.X.data(data.X.data > 180) = data.X.data(data.X.data > 180) - 360;

ncid = netcdf.open(hycom.Latitude, 'NOWRITE');
try
    varid = netcdf.inqVarID(ncid, suffix.Latitude{1});
catch
    varid = netcdf.inqVarID(ncid, suffix.Latitude{2});
end

data.Y.data = netcdf.getVar(ncid, varid, 'double');

netcdf.close(ncid)

% If the spatial data are vectors, turn them in to matrices here.
if isvector(data.X.data) && isvector(data.Y.data)
    [data.X.data, data.Y.data] = meshgrid(data.X.data, data.Y.data);
    % Orient the arrays as required for the calls to the remote server.
    data.X.data = data.X.data';
    data.Y.data = data.Y.data';
end

% Create indices of the size of the arrays.
data.X.idx = repmat(1:size(data.X.data, 1), [size(data.X.data, 2), 1])';
data.Y.idx = repmat(1:size(data.Y.data, 2), [size(data.Y.data, 1), 1]);
%data.Y.idx = 1:size(data.Y.data, 2) - 1;

% Find the indices which cover the model domain then find the extremes to
% request only a subset from the OPeNDAP server.
idx = data.X.data > extents(1) & data.X.data < extents(2) & data.Y.data > extents(3) & data.Y.data < extents(4);
xrange = [min(data.X.idx(idx)), max(data.X.idx(idx))];
yrange = [min(data.Y.idx(idx)), max(data.Y.idx(idx))];

data.lon = data.X.data(xrange(1):xrange(2), yrange(1):yrange(2));
data.lat = data.Y.data(xrange(1):xrange(2), yrange(1):yrange(2));

% Clear out the full lon/lat arrays.
% data = rmfield(data, {'X', 'Y'});

% Now get the variables for which we've been asked.
fields = varlist;

for aa = 1:length(fields)
    % Store the downloaded data in a struct. Assume the spatial
    % data is identical to that in data.lon and data.lat.
    data.(fields{aa}).data = [];

    ncid = netcdf.open(hycom.(fields{aa}));
    varid = netcdf.inqVarID(ncid, suffix.(fields{aa}){name_index});

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
        url = get_url(times(tt), threehourly);

        if ~strcmpi(oldurl, url) || tt == 1
            if times(tt) < greg2mjulian(2008, 09, 19, 0, 0, 0)
                hycom.(fields{aa}) = [url, suffix.(fields{aa}){1}];
            elseif times(tt) >= greg2mjulian(2008, 09, 19, 0, 0, 0)
                hycom.(fields{aa}) = [url, suffix.(fields{aa}){2}];
            end
            % Close any existing open connections and reopen with
            % the new URL.
            netcdf.close(ncid)
            ncid = netcdf.open(hycom.(fields{aa}));
            try
                varid = netcdf.inqVarID(ncid, suffix.(fields{aa}){1});
            catch
                varid = netcdf.inqVarID(ncid, suffix.(fields{aa}){2});
            end

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
    netcdf.close(ncid);
end

if ftbverbose
    fprintf('end   : %s\n', subname)
end

function url = get_url(time, threehourly)
% Child function to find the relevant URL to use for a given time step.
%
% url = get_url(time, threehourly);
%
% INPUT:
%   time - Modified Julian Day
%   threehourly - additional string to append for optional three hourly
%   outputs for the 1992/10/02 to 2008/09/18 period. Leave empty for daily
%   data.
%
% OUTPUT:
%   url - string of the approprate URL for the date supplied in time.
%

[t1, t2, t3, t4, t5, t6] = datevec(datestr(now));

if time < greg2mjulian(1992, 10, 2, 0, 0, 0)
    error('No HYCOM data available prior to 1992-10-02. Select a start date from 1992-10-02 onwards.')
elseif time >= greg2mjulian(1992, 10, 2, 0, 0, 0) && time < greg2mjulian(1995, 7, 31, 0, 0, 0)
    warning('Using the HYCOM Global Reanalysis data for dates up to 2008/09/16, thereafter the Global Analysis.')
    url = sprintf('http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.0%s?', threehourly);
elseif time >= greg2mjulian(1995, 7, 31, 0, 0, 0) && time < greg2mjulian(2008, 09, 19, 0, 0, 0)
    warning('Using the HYCOM Global Reanalysis data for dates up to 2008/09/16, thereafter the Global Analysis.')
    url = sprintf('http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1%s?', threehourly);
elseif time >= greg2mjulian(2008, 9, 19, 0, 0, 0) && time < greg2mjulian(2009, 5, 7, 0, 0, 0)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.6?';
elseif time >= greg2mjulian(2009, 5, 7, 0, 0, 0) && time < greg2mjulian(2011, 1, 3, 0, 0, 0)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.8?';
elseif time >= greg2mjulian(2011, 1, 3, 0, 0, 0) && time < greg2mjulian(2013, 8, 21, 0, 0, 0)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9?';
elseif time >= greg2mjulian(2013, 8, 21, 0, 0, 0) && time < greg2mjulian(2014, 4, 21, 0, 0, 0):
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.0?';
elseif time >= greg2mjulian(2014, 4, 21, 0, 0, 0) && time < greg2mjulian(2016, 4, 18, 0, 0, 0):
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1?';
elseif time >= greg2mjulian(2016, 4, 18, 0, 0, 0) && time <= greg2mjulian(t1, t2, t3, t4, t5, t6):
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.2?';
elseif time > greg2mjulian(t1, t2, t3, t4, t5, t6)
    error('Given date is in the future.')
else
    error('Date is outside of the known spacetime continuum. See help TARDIS.')
end

