function data = get_HYCOM_forcing(Mobj, modelTime)
% Get mean flow, temperature and salinity data from HYCOM model outputs
% through their OPeNDAP server. 
%
% data = get_HYCOM_forcing(Mobj, modelTime)
%
% DESCRIPTION:
%   Using OPeNDAP, extract the necessary parameters to create an FVCOM
%   forcing file. Requires the OPeNDAP toolbox (see below for where to get
%   it) for versions of MATLAB older than 2012a (v7.14).
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
%     - u flow component
%     - v flow component
%
% REQUIRES:
%   The OPeNDAP toolbox:
%       http://www.opendap.org/pub/contributed/source/ml-toolbox/
%
% Author(s)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2013-01-31 First version.
%   2013-08-19 Make some more progress in getting this working. Mainly
%   change the way the dates are handled to form the relevant URL for
%   downloading the data.
%
%==========================================================================

subname = 'get_HYCOM_forcing';

global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

% For checking whether to use the third-party OPeNDAP toolbox or native
% MATLAB tools. OPeNDAP support is included in MATLAB version 7.14 onwards.
v714date = datenum(2012, 3, 1);
currdate = ver('MATLAB');
currdate = datenum(currdate.Date);

% Get the extent of the model domain (in spherical)
if ~Mobj.have_lonlat
    error('Need spherical coordinates to extract the forcing data')
else
    % Add a buffer of one grid cell in latitude and two in longitude to
    % make sure the model domain is fully covered by the extracted data.
    [dx, dy] = deal(1/12, 1/12); % HYCOM resolution in degrees
    extents = [min(Mobj.lon(:))-(2*dx), max(Mobj.lon(:))+(2*dx), min(Mobj.lat(:))-dy, max(Mobj.lat(:))+dy];
end

if modelTime(end) - modelTime(1) > 365
    error('Can''t (yet) process more than a year at a time.')
end

[yearStart, monStart, dayStart, hStart, mStart, sStart] = mjulian2greg(modelTime(1));
[yearEnd, monEnd, dayEnd, hEnd, mEnd, sEnd] = mjulian2greg(modelTime(end));

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

% Get the days of year for the start and end dates
indStart = floor(datenum([yearStart, monStart, dayStart, hStart, mStart, sStart])) - datenum([yearStart, 1, 1]);
indEnd = ceil(datenum([yearEnd, monEnd, dayEnd, hEnd, mEnd, sEnd])) - datenum([yearEnd, 1, 1]);
% Create the necessary string for the time indices.
tInd = sprintf('[%i:1:%i]', indStart, indEnd);

% Set up a struct of the HYCOM data sets in which we're interested.
hycom.Longitude     = [url, '?Longitude'];
hycom.Latitude      = [url, '?Latitude'];
hycom.temperature   = [url, '?temperature'];
hycom.salinity      = [url, '?salinity'];
% hycom.density       = [url, '?density'];
hycom.u             = [url, '?u']; % mean flow
hycom.v             = [url, '?v']; % mean flow
hycom.MT            = [url, '?MT']; % time
% hycom.X             = [url, '?X']; % causes MATLAB to crash...
% hycom.Y             = [url, '?Y']; % causes MATLAB to crash...

% % The HYCOM OPeNDAP server is spectacularly slow (at the moment?). As such,
% % extracting even the entire lat/long arrays is painfully slow. Instead of
% % directly interrogating the server, we'll assume a uniformly spaced 1/12th
% % degree grid and use that to find the indices. We'll give ourself a bit of
% % wiggle-room to accommodate the inevitable inaccuracy in this method.
% [grid.X, grid.Y] = meshgrid(0:1/12:360, -90:1/12:90);
% grid.ids = inbox(extents, grid.X, grid.Y);
% 
% if extents(1) < 0 && extents(2) < 0
%     % This is OK, we can just shunt the values by 360.
%     extents(1) = extents(1) + 360;
%     extents(2) = extents(2) + 360;
%     index_lon = find(data_lon.lon > extents(1) & data_lon.lon < extents(2));
% elseif extents(1) < 0 && extents(2) > 0
%     % This is the tricky one. We'll do two passes to extract the
%     % western chunk first (extents(1)+360 to 360), then the eastern
%     % chunk (0-extents(2))
%     index_lon{1} = find(data_lon.lon >= extents(1) + 360);
%     index_lon{2} = find(data_lon.lon <= extents(2));
% else
%     % Dead easy, we're in the eastern hemisphere, so nothing too
%     % strenuous here
%     index_lon = find(data_lon.lon > extents(1) & data_lon.lon < extents(2));
% end
% 
% % Latitude is much more straightforward
% data_lat = loaddap([ncep.(fields{aa}),'?lat']);
% index_lat = find(data_lat.lat > extents(3) & data_lat.lat < extents(4));

% Now get the latitude and longitude to calculate the indices representing
% the model domain.

% Depending on the MATLAB version, use either the third-party OPeNDAP
% toolbox or the native MATLAB tools.

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
    
    % Make the longitude values in the range 0-360.
    data.X.data = mod(data.X.data, 360);

    ncid = netcdf.open(hycom.Latitude, 'NOWRITE');
    varid = netcdf.inqVarID(ncid, 'Latitude');

    data.Y.data = netcdf.getVar(ncid, varid, 'double');

    % Due to the dual polar nature of the HYCOM model domain, and since I'm
    % lazy, I'm going to throw away the data above the rectilinear part of
    % the model (4500x2172 rather than 4500x3298).
    data.X.data = data.X.data(:, 
    data.X.idx = 1:size(data.X.data, 1) - 1;
    data.Y.idx = 1:size(data.X.data, 1) - 1;
    
    % Find the indices which fit within the model domain.
    data.X.idx = data.X.idx
    
else
    data.X.idx = loaddap(hycom.X);
    data.Y.idx = loaddap(hycom.Y);
    xIdx = length(data.X.idx.X) - 1;
    yIdx = length(data.Y.idx.Y) - 1;
    data.lon.all = loaddap([hycom.lon, sprintf('[%i:1:%i]', 0, 0), sprintf('[%i:1:%i]', 0, xIdx)]);
    data.lat.all = loaddap([hycom.lat, sprintf('[%i:1:%i]', 0, yIdx), sprintf('[%i:1:%i]', 0, 0)]);
end

fields = fieldnames(hycom);

for aa = 1:length(fields)   
    % Store the downloaded data in a struct with associated spatial and
    % temporal data.
    data.(fields{aa}).data = [];
    data.(fields{aa}).time = [];
    data.(fields{aa}).lat = [];
    data.(fields{aa}).lon = [];
    data_attributes.(fields{aa}) = [];

    % Get attributes from which to calculate length of various fields.
    if datenum(currdate) > v714date
        
        ncid = netcdf.open(hycom.(fields{aa}));

        % If you don't know what it contains, start by using the
        % 'netcdf.inq' operation:
        %[numdims,numvars,numglobalatts,unlimdimid] = netcdf.inq(ncid);

        varid = netcdf.inqVarID(ncid, fields{aa});
        
        data.(fields{aa}).data = netcdf.getVar(ncid, varid, 'double');
        
        
       
        data_attributes.(fields{aa}) = loaddap('-A', [hycom.(fields{aa})]);
        % Get the data time and convert to Modified Julian Day.
        data_time = loaddap(hycom.time);
    else
        data_attributes.(fields{aa}) = loaddap('-A', [hycom.(fields{aa})]);
        % Get the data time and convert to Modified Julian Day.
        data_time = loaddap(hycom.time);
    end

    timevec = datevec((data_time.MT)/24+365);
    data.time = greg2mjulian(timevec(:,1), timevec(:,2), timevec(:,3), ...
        timevec(:,4), timevec(:,5), timevec(:,6));
    % Clip the time to the given range
    data_time_mask = data.time >= modelTime(1) & data.time <= modelTime(end);
    data_time_idx = 1:size(data.time,1);
    data_time_idx = data_time_idx(data_time_mask);
    data.time = data.time(data_time_mask);
end
