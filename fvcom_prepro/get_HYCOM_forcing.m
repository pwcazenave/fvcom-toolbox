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

% Get the extent of the model domain (in spherical)
if ~Mobj.have_lonlat
    error('Need spherical coordinates to extract the forcing data')
else
    % Add a buffer of two grid cells in latitude and two in longitude to
    % make sure the model domain is fully covered by the extracted data.
    [dx, dy] = deal(1/12, 1/12); % HYCOM resolution in degrees
    % West, east, south, north
    extents = [min(Mobj.lon(:))-(2*dx), max(Mobj.lon(:))+(2*dx), ...
        min(Mobj.lat(:))-(2*dy), max(Mobj.lat(:))+(2*dy)];
end

if modelTime(end) - modelTime(1) > 365
    error('Can''t (yet) process more than a year at a time.')
end

[yearStart, ~, ~, ~, ~, ~] = mjulian2greg(modelTime(1));
[yearEnd, ~, ~, ~, ~, ~] = mjulian2greg(modelTime(end));

if yearEnd ~= yearStart
    error('Can''t (yet) process across a year boundary.')
else
    year = yearEnd;
end

% Use Modified Julian Days for the time checks
[t1, t2, t3, t4, t5, t6] = datevec(date);
if modelTime(1) < greg2mjulian(2008, 09, 16, 0, 0, 0)
    error('Not using the legacy HYCOM model output. Select a start date from September 2008 onwards.')
elseif modelTime(1) >= greg2mjulian(2008, 9, 16, 0, 0, 0) && modelTime(1) < greg2mjulian(2009, 5, 6, 0, 0, 0)
    url = ['http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.6', num2str(year)];
elseif modelTime(1) >= greg2mjulian(2009, 5, 6, 0, 0, 0) && modelTime(1) < greg2mjulian(2011, 1, 2, 0, 0, 0)
    url = ['http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.8', num2str(year)];
elseif modelTime(1) >= greg2mjulian(2011, 1, 2, 0, 0, 0) && modelTime(1) <= greg2mjulian(t1, t2, t3, 0, 0, 0)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9';
elseif modelTime(1) > greg2mjulian(t1, t2, t3,  t4, t5, t6)
    error('Date is in the future?')
else
    error('Date is outside of known spacetime?')
end
% REMOVE ME!!!
url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9';
% REMOVE ME!!!

% Set up a struct of the HYCOM data sets in which we're interested.
hycom.Longitude     = [url, '?Longitude'];
hycom.Latitude      = [url, '?Latitude'];
hycom.temperature   = [url, '?temperature'];
hycom.salinity      = [url, '?salinity'];
% hycom.density       = [url, '?density']; % don't need for now
hycom.u             = [url, '?u']; % mean flow
hycom.v             = [url, '?v']; % mean flow
hycom.Depth         = [url, '?Depth']; % water depth
% hycom.MT            = [url, '?MT']; % time - don't need to get automatically
% hycom.X             = [url, '?X']; % causes MATLAB to crash...
% hycom.Y             = [url, '?Y']; % causes MATLAB to crash...

data.X.data = [];
data.Y.data = [];
data.time.data = [];

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

    % Now find the time indices. HYCOM stores its times as days since
    % 1900-12-31 00:00:00. Naturally this is completely different from
    % everything else ever.
    ncid = netcdf.open(hycom.MT, 'NOWRITE');
    varid = netcdf.inqVarID(ncid, 'MT');
    data.MT.data = netcdf.getVar(ncid, varid, 'double');
    % Convert to MATLAB days and then to Modified Julian Days.
    t = datevec(data.MT.data + datenum(1900, 12, 31, 0, 0, 0));
    tj = greg2mjulian(t(:,1), t(:,2), t(:,3), t(:,4), t(:,5), t(:,6));
    % Find the time indices which cover the model period.
    trange = [find(tj >= modelTime(1), 1, 'first'), find(tj <= modelTime(2), 1, 'last')];

    data.time = tj; % Export the Modified Julian Day time series.

    % Clear out the full lon/lat arrays.
    data = rmfield(data, {'X', 'Y'});

else
    % I haven't tested this at all.

    data.X.idx = loaddap(hycom.X);
    data.Y.idx = loaddap(hycom.Y);
    xIdx = length(data.X.idx.X) - 1;
    yIdx = length(data.Y.idx.Y) - 1;
    data.lon.all = loaddap([hycom.lon, sprintf('[%i:1:%i]', 0, 0), sprintf('[%i:1:%i]', 0, xIdx)]);
    data.lat.all = loaddap([hycom.lat, sprintf('[%i:1:%i]', 0, yIdx), sprintf('[%i:1:%i]', 0, 0)]);
end

netcdf.close(ncid);

fields = fieldnames(hycom);

for aa = 1:length(fields)

    if strcmpi(fields{aa}, 'Longitude') || strcmpi(fields{aa}, 'Latitude')
        % Skip these as we've already got the spatial information as
        % data.lon and data.lat.
        continue
    else

        % Store the downloaded data in a struct. Assume the spatial data is
        % identical to that in data.lon and data.lat.
        data.(fields{aa}).data = [];
        data.(fields{aa}).time = [];
        data_attributes.(fields{aa}) = [];

        % Get attributes from which to calculate length of various fields.
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
            % extract a day's data at a time and stick them together here.
            % If nothing else, this at least means I can give some
            % indication of progress rather than just wait for something to
            % magically happen.
            ts = trange(1):trange(2);
            nt = length(ts);
            nx = (xrange(2) - xrange(1)) + 1;
            ny = (yrange(2) - yrange(1)) + 1;
            nz = 33; % hard code vertical levels

            data.(fields{aa}).data = nan(nx, ny, nz, nt);
            for tt = 1:nt
                if ftbverbose
                    fprintf('%s: time %i of %i... ', fields{aa}, tt, nt)
                end
                start = [xrange(1), yrange(1), 1, ts(tt)] - 1;
                count = [(xrange(2) - xrange(1)) + 1, (yrange(2) - yrange(1)) + 1, nz, 1];

                data.(fields{aa}).data(:, :, :, tt) = netcdf.getVar(ncid, varid, start, count, 'double');
                if ftbverbose
                    fprintf('done.\n')
                end
            end

        else
            % Third-party toolbox version (unimplemented).

            data_attributes.(fields{aa}) = loaddap('-A', [hycom.(fields{aa})]);
            % Get the data time and convert to Modified Julian Day.
            data_time = loaddap(hycom.time);
        end
    end
end

netcdf.close(ncid);