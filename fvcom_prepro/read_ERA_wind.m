function era = read_ERA_wind(year, datadir, varlist)
% Reads in ERA Interim files and outputs a struct containing the requested
% variables for a given year. 
% 
% ERA = read_ERA_wind(YEAR, DATADIR, VARLIST)
% 
% DESCRIPTION:
%   For the given YEAR, find all the ERA Interim NetCDF files and aggregate
%   them into a single MATLAB struct, which contains the variables
%   specified in VARLIST. In addition to the specified variables, time and
%   longitude and latitude will also be extracted. 
% 
% INPUT:
%   YEAR - year to extract
%   DATADIR - path to the directory which contains the ERA NetCDF files
%   VARLIST - list of the particular variables to extract from the NetCDF
%   files.
% 
% OUTPUT:
%   era - struct containing the variables specified in VARLIST.
% 
% Author(s)
%   Pierre Cazenave (Plymouth Marine Laboratory)
% 
% Revision history:
%   2012-10-19 First version based loosely on read_NCEP_wind.m in the
%   fvcom-toolbox.
% 
%==========================================================================

% year = 2006;
% % datadir = '/data/modellers/to_archive/momm-ERA40-interim/';
% datadir = '/users/modellers/pica/Data/ECMWF/2006';
% varlist = {'u10', 'v10'};

warning off

if nargin ~= 3
    error('Incorrect number of arguments')
end

subname = 'read_ERA_wind';

global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

if exist(datadir, 'dir') ~= 7
   error(['file: ' datadir ' does not exist']);
end

%--------------------------------------------------------------------------
% Open ERA Interim data and check for time range
%--------------------------------------------------------------------------

% Get the time. 
ncERA = netcdf.open(fullfile(datadir, [num2str(year), '_U10V10.nc']), 'NOWRITE');
time_varid = netcdf.inqVarID(ncERA, 'time');
eratimehours = netcdf.getVar(ncERA, time_varid);

% ERA Interim times are stored as hours since 1900-01-01 00:00:0.0.
% MATLAB's dates days since 0000/00/00 00:00:00. 
eratimehours = datevec((double(eratimehours)/24) + datenum('1900-01-01 00:00:0.0'));
% Convert the ERA times to Modified Julian Date.
era.time = greg2mjulian(eratimehours(:,1), eratimehours(:,2),...
    eratimehours(:,3), eratimehours(:,4), eratimehours(:,5),...
    eratimehours(:,6));

if(ftbverbose);
    fprintf('beg time of ERA Interim data %04i/%02i/%02i %02i:%02i:%02i\n', eratimehours(1,:));
    fprintf('end time of ERA Interim data %04i/%02i/%02i %02i:%02i:%02i\n', eratimehours(end,:));
end

% Get the geographical information from the ERA Interim data. Again, use
% the U10 file only (we're assuming they're both global).
lat_varid = netcdf.inqVarID(ncERA, 'latitude');
lon_varid = netcdf.inqVarID(ncERA, 'longitude');
eralatvector = netcdf.getVar(ncERA, lat_varid);
eralonvector = netcdf.getVar(ncERA, lon_varid);
[era.lon, era.lat] = meshgrid(eralonvector, eralatvector);

% Find the necessary variables
for var=1:numel(varlist)

    getVar = varlist{var};
    varid_ERA = netcdf.inqVarID(ncERA, getVar);

    % Get the data
    data = netcdf.getVar(ncERA, varid_ERA, 'single');
    
    if strcmpi(getVar, 'u10') || strcmpi(getVar, 'v10')
        % The ERA Interim wind component data are packed as integers. The
        % following equation describes how to unpack them:
        %     unpacked value = add_offset + ((packed value)*scale_factor)
        % (from
        % http://www.ecmwf.int/products/data/archive/data_faq.html#netcdfintegers).
        % ERA wind scale_factor and add_offset are doubles (the NCEP ones
        % are singles).
        scale_factor = netcdf.getAtt(ncERA,varid_ERA,'scale_factor','double');
        add_offset = netcdf.getAtt(ncERA,varid_ERA,'add_offset','double');

        % Unpack the values. In general, the data for U10 and V10 should be
        % doubles for griddata to work. Fix the order of the dimensions to
        % match the coordinates in eralon and eralat.
        era.(getVar) = permute(double(add_offset + (data.*scale_factor)), [2,1,3]);
    else
        % We're assuming they're not packed and so we just return the data
        % as is (but as doubles).
        era.(getVar) = double(data);
    end
end

netcdf.close(ncERA)

if ftbverbose;
    fprintf(['end   : ' subname '\n'])
end