function ncep = read_NCEP_wind(ncep_u10_file, ncep_v10_file)
% Reads in two NCEP wind vector files (U and V) and outputs four arrays of
% longitude, latitude, u10 and v10 velocity components.
% 
% function read_NCEP_wind()
% 
% DESCRIPTION:
%   Read a pair of NCEP NetCDF files (U10 and V10 vectors) and output to
%   four arrays of longitude, latitude, u10 and v10.
% 
% INPUT:
%   NCEP NetCDF U10 filename (and path)
%   NCEP NetCDF V10 filename (and path)
% 
% OUTPUT:
%   ncep - struct with the time, latitude, longitude, u10 and v10 arrays in
%   it. Time is in Modified Julian Days.
% 
% Author(s)
%   Pierre Cazenave (Plymouth Marine Laboratory)
% 
% Revision history:
%   2012-10-16 First version based on parts of ncep2fvcom_U10V10.m in the
%   fvcom-toolbox.
% 
%==========================================================================

warning off

if nargin ~= 2
    error('Incorrect number of arguments')
end

subname = 'read_NCEP_wind';

global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

if exist(ncep_u10_file, 'file') ~= 2
   error(['file: ' ncep_u10_file ' does not exist']);
end
if exist(ncep_v10_file, 'file') ~= 2
    error(['file: ' ncep_v10_file ' does not exist']);
end

%--------------------------------------------------------------------------
% Open NCEP data and check for time range
%--------------------------------------------------------------------------

% Get the year from the NCEP file name
ncep_u10_year = get_NCEP_year(ncep_u10_file);
ncep_v10_year = get_NCEP_year(ncep_v10_file);
% Check both files are from the same year
if (ncep_u10_year - ncep_v10_year) ~= 0
    error('Input U and V wind data files are from different years')
end

% Get the time. It's stored relative to the 1st January for the given year.
% We're assuming that since both files are for the same year, we don't have
% to pull the 'time' variable out from both (they should be identical).
nc_u10 = netcdf.open(ncep_u10_file, 'NOWRITE');
nc_v10 = netcdf.open(ncep_v10_file, 'NOWRITE');
time_varid = netcdf.inqVarID(nc_u10, 'time');
nceptimehours = netcdf.getVar(nc_u10, time_varid);

% NCEP dates are relative to 0001/01/01 00:00:00 and stored in hours.
% MATLAB's dates are relative to 0000/00/00 00:00:00 and stored in days.
% Need to add a year and a day to the NCEP time when converting.
nceptimedays = datevec((nceptimehours/24) + datenum(1, 0, -1));
ncep.time = greg2mjulian(nceptimedays(:,1), nceptimedays(:,2),...
    nceptimedays(:,3), nceptimedays(:,4), nceptimedays(:,5),...
    nceptimedays(:,6));

if(ftbverbose);
    fprintf('beg time of NCEP data %04i/%02i/%02i %02i:%02i:%02i\n',nceptimedays(1,:));
    fprintf('end time of NCEP data %04i/%02i/%02i %02i:%02i:%02i\n',nceptimedays(end,:));
end

% Get the geographical information from the NCEP data. Again, use the U10
% file only (we're assuming they're both global).
lat_varid = netcdf.inqVarID(nc_u10, 'lat');
lon_varid = netcdf.inqVarID(nc_u10, 'lon');
nceplatvector = netcdf.getVar(nc_u10, lat_varid);
nceplonvector = netcdf.getVar(nc_u10, lon_varid);

[ncep.lon, ncep.lat] = meshgrid(nceplonvector, nceplatvector);

% Find the necessary variables
u10_varid_NCEP = netcdf.inqVarID(nc_u10, 'uwnd');
v10_varid_NCEP = netcdf.inqVarID(nc_v10, 'vwnd');

% Get the U10 and V10 data
U10 = netcdf.getVar(nc_u10, u10_varid_NCEP, 'single');
V10 = netcdf.getVar(nc_v10, v10_varid_NCEP, 'single');

% The NCEP data are packed as integers. The following equation describes
% how to unpack them:
%     unpacked value = add_offset + ( (packed value) * scale_factor )
% (from http://www.esrl.noaa.gov/psd/data/gridded/faq.html#2).
% Keep them as singles for now to avoid horrible rounding errors.
scale_factor = netcdf.getAtt(nc_u10,u10_varid_NCEP,'scale_factor','single');
add_offset = netcdf.getAtt(nc_u10,u10_varid_NCEP,'add_offset','single');

% Unpack the values. U10 and V10 must be doubles for griddata to work. Fix
% the order of the dimensions to match the coordinates in nceplon and
% nceplat. 
ncep.uwnd = permute(double(add_offset + (U10.*scale_factor)), [2,1,3]);
ncep.vwnd = permute(double(add_offset + (V10.*scale_factor)), [2,1,3]);

netcdf.close(nc_u10)
netcdf.close(nc_v10)
