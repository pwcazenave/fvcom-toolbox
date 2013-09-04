function data = get_HYCOM_forcing(Mobj, modelTime)
% Get mean flow, temperature and salinity data from HYCOM model outputs
% through their OPeNDAP server. 
%
% data = get_HYCOM_forcing(Mobj, modelTime)
%
% DESCRIPTION:
%   Using OPeNDAP, extract the necessary parameters to create an FVCOM
%   forcing file.
%
% INPUT: 
%   Mobj - MATLAB mesh object
%   modelTime - Modified Julian Date start and end times
%
% OUTPUT:
%   data - struct of the data necessary to force FVCOM. These can be
%   interpolated onto an unstructured grid in Mobj using
%   grid2fvcom.m.
%
% The parameters which are obtained from the HYCOM data are:
%     - temperature
%     - salinity
%     - u mean flow component
%     - v mean flow component
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
%   moment, I'm ignoring the old toolbox for reading the data from the
%   HYCOM OPeNDAP server.
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
suffix.X            = '?X';
suffix.Y            = '?Y';

% Set up a struct of the HYCOM data sets in which we're interested.
hycom.MT            = [url, suffix.MT];             % time [1D]
hycom.Longitude     = [url, suffix.Longitude];      % [2D]
hycom.Latitude      = [url, suffix.Latitude];       % [2D]
hycom.Depth         = [url, suffix.Depth];          % water depth [2D]
hycom.temperature   = [url, suffix.temperature];    % [4D]
hycom.salinity      = [url, suffix.salinity];       % [4D]
hycom.u             = [url, suffix.u];              % mean flow % [4D]
hycom.v             = [url, suffix.v];              % mean flow % [4D]
% hycom.density       = [url, suffix.density];        % don't need for now
% hycom.X             = [url, suffix.X];              % crashes MATLAB...
% hycom.Y             = [url, suffix.Y];              % crashes MATLAB...

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

            % Convert to MATLAB days and then to Modified Julian Days.
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

fields = fieldnames(hycom);

for aa = 1:length(fields)

    switch fields{aa}

        case {'Longitude', 'Latitude', 'MT', 'Depth'}
            % Skip these as we've already got the information at the
            % beginning of this function. These are largely time
            % independent data, so no need to get them here as they don't
            % have 4 dimensions, whereas the code after the otherwise
            % assumed 4D data.
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
                % We probably want all depths but only a subset in time and
                % space.

                % Since the HYCOM OPeNDAP server is so spectacularly slow,
                % extract a day's data at a time and stick them together
                % here. If nothing else, this at least means I can give
                % some indication of progress rather than just wait for
                % something to magically happen.
                nx = (xrange(2) - xrange(1)) + 1;
                ny = (yrange(2) - yrange(1)) + 1;

                % Preallocate the output so we don't append to an array
                % (which is slow).
                data.(fields{aa}).data = nan(nx, ny, nz, nt);

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

                    % netCDF starts at zero, hence -1.
                    start = [xrange(1), yrange(1), 1, ts] - 1;
                    count = [nx, ny, nz, 1];

                    data.(fields{aa}).data(:, :, :, tt) = netcdf.getVar(ncid, varid, start, count, 'double');

                    if ftbverbose; fprintf('done.\n'); end
                end
            else
                % Third-party toolbox version (unimplemented).

                data_attributes.(fields{aa}) = loaddap('-A', [hycom.(fields{aa})]);
                % Get the data time and convert to Modified Julian Day.
                data.time = loaddap(hycom.time);
            end
    end
end

netcdf.close(ncid);


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
