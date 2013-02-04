function data = get_HYCOM_forcing(Mobj, modelTime)
% Get mean flow, temperature and salinity data from HYCOM model outputs
% through their OPeNDAP server. 
% 
% data = get_HYCOM_forcing(Mobj, modelTime)
% 
% DESCRIPTION:
%   Using OPeNDAP, extract the necessary parameters to create an FVCOM
%   forcing file. Requires the OPeNDAP toolbox (see below for where to get
%   it).
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
%     - sea surface height (ssh)
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
% 
%==========================================================================

subname = 'get_HYCOM_forcing';

global ftbverbose;
if(ftbverbose);
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end

% Get the extent of the model domain (in spherical)
if ~Mobj.have_lonlat
    error('Need spherical coordinates to extract the forcing data')
else
    % Add a buffer of one grid cell in latitude and two in longitude to
    % make sure the model domain is fully covered by the extracted data.
    [dx, dy] = deal(2.5, 2.5); % NCEP resolution in degrees
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

if year < 2008 && monStart < 9
    error('Not using the legacy HYCOM model output. Select a start date from September 2008 onwards.')
elseif (year >= 2008 && monStart >= 9) && (year < 2009 && monEnd <= 5)
    url = ['http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.6/', num2str(year)];
elseif (year >= 2009 && monStart >= 5) && (year < 2011 && monEnd <= 1)
    url = ['http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.8/', num2str(year)];
elseif (year >= 2011 && monStart >= 1)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9';
else
    error('Date is bef?')
end

% Get the days of year for the start and end dates
indStart = floor(datenum([yearStart, monStart, dayStart, hStart, mStart, sStart])) - datenum([yearStart, 1, 1]);
indEnd = ceil(datenum([yearEnd, monEnd, dayEnd, hEnd, mEnd, sEnd])) - datenum([yearEnd, 1, 1]);
% Create the necessary string for the time indices.
tInd = sprintf('[%i:1:%i]', indStart, indEnd);

% Set up a struct of the HYCOM data sets in which we're interested.
hycom.temp  = [url, '?temperature'];
hycom.salt  = [url, '?salinity'];
hycom.u     = [url, '?u'];
hycom.v     = [url, '?v'];
hycom.ssh   = [url, '?ssh'];
hycom.time  = [url, '?MT'];
hycom.X     = [url, '?X'];
hycom.Y     = [url, '?Y'];
hycom.lon   = [url, '?Longitude'];
hycom.lat   = [url, '?Latitude'];

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
data.X.idx = loaddap(hycom.X);
data.Y.idx = loaddap(hycom.Y);
xIdx = length(data.X.idx.X) - 1;
yIdx = length(data.Y.idx.Y) - 1;
data.lon.all = loaddap([hycom.lon, sprintf('[%i:1:%i]', 0, 0), sprintf('[%i:1:%i]', 0, xIdx)]);
data.lat.all = loaddap([hycom.lat, sprintf('[%i:1:%i]', 0, yIdx), sprintf('[%i:1:%i]', 0, 0)]);

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
    data_attributes.(fields{aa}) = loaddap('-A', [hycom.(fields{aa})]);
    % Get the data time and convert to Modified Julian Day.
    data_time = loaddap(hycom.time);

    timevec = datevec((data_time.time)/24+365);
    data.time = greg2mjulian(timevec(:,1), timevec(:,2), timevec(:,3), ...
        timevec(:,4), timevec(:,5), timevec(:,6));
    % Clip the time to the given range
    data_time_mask = data.time >= modelTime(1) & data.time <= modelTime(end);
    data_time_idx = 1:size(data.time,1);
    data_time_idx = data_time_idx(data_time_mask);
    data.time = data.time(data_time_mask);

