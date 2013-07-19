function data = get_MetUM_forcing(Mobj, modelTime, credentials)
% Get the required parameters from the Met Office Unified Model (TM)
% (hereafter MetUM) to use in FVCOM surface forcing.
%
% data = get_MetUM_forcing(Mobj, modelTime, credentials)
%
% DESCRIPTION:
%   Using FTP access, extract the necessary parameters to create an FVCOM
%   forcing file. Requires the air_sea toolbox (see below for where to get
%   it). Data are sampled four times daily.
%
% INPUT:
%   Mobj - MATLAB mesh object
%   modelTime - Modified Julian Date start and end times
%   credentials - struct with fields username and password to access the
%   FTP server.
%
% OUTPUT:
%   data - struct of the data necessary to force FVCOM. These can be
%   interpolated onto an unstructured grid in Mobj using grid2fvcom.m.
%
% The required parameters which can be obtained are:
%     - surface_net_downward_shortwave_flux (W m-2)
%     - surface_downwelling_shortwave_flux_in_air (W m-2)
%     - surface_net_downward_longwave_flux (W m-2)
%     - surface_downwelling_longwave_flux_in_air (W m-2)
%     - surface_upward_sensible_heat_flux (W m-2)
%     - eastward_wind / x_wind (m s-1)
%     - northward_wind / y_wind (m s-1)
%     - surface_upward_latent_heat_flux  (W m-2)
%     - air_temperature (K)
%     - relative_humidity (%)
%     - precipitation_flux (kg m-2 s-1)
%     - air_pressure_at_sea_level (Pa)
%
% In addition to these, the momentum flux is calculated from wind data.
% Precipitation is converted from kg/m^2/s to m/s. Evaporation is
% calculated from the mean daily latent heat net flux (lhtfl) at the
% surface.
%
% REQUIRES:
%   The air_sea toolbox:
%       http://woodshole.er.usgs.gov/operations/sea-mat/air_sea-html/index.html
%
% Author(s)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2013-05-07 First version.
%
%==========================================================================
subname = 'get_MetUM_forcing';

global ftbverbose;
if ftbverbose;
    fprintf('\nbegin : %s \n', subname)
end

% Get the extent of the model domain (in spherical)
if ~Mobj.have_lonlat
    error('Need spherical coordinates to extract the forcing data')
else
    % Add a 1 degree buffer to make sure the model domain is fully covered
    % by the extracted data.
    [dx, dy] = deal(1, 1);
    extents = [min(Mobj.lon(:))-(2*dx), max(Mobj.lon(:))+(2*dx), min(Mobj.lat(:))-dy, max(Mobj.lat(:))+dy];
end

nt = modelTime(end) - modelTime(1);
if nt > 365
    error('Can''t (yet) process more than a year at a time.')
end

[yearStart, monthStart, dayStart] = mjulian2greg(modelTime(1));
[yearEnd, monthEnd, dayEnd] = mjulian2greg(modelTime(end));
t = modelTime(1):1/4:modelTime(end);

if yearEnd ~= yearStart
    error('Can''t (yet) process across a year boundary.')
end

if yearStart < 2006 || yearEnd > 2012
    error('The MetUM repository does not contain data earlier than 2006 and later than 2012')
end

% For the pre-2006 data, we need to download several files with unique
% names. The names are based on the STASH numbers and the date:
%   naamYYYYMMDDHH_STASH#_00.pp
% The numbers we're interested in are stored in stash.
stash = [2, 3, 407, 408, 409, 4222, 9229, 16004];
vars = {'uwnd', 'uwnd', 'vwnd', 'vwnd', 'slp_rho', 'slp_theta', ...
    'surface_air_pressure', 'air_sw', 'air_lw', ...
    'rhum', 'prate', 'temp_model', 'temp_press'};

ns = length(stash);

% From where will we be downloading the data?
site = 'ftp.ceda.ac.uk';
basePath = 'badc/ukmo-um/data/nae/';

% Open a remote connection to the FTP site
remote = ftp(site, credentials(1), credentials(2));

% Depending on the year we're extracting, we need to append different
% directories to get the data. 
for i = 1:nt * 4 % four files per day (at 0000, 0600, 1200 and 1800).

    [year, month, day, hour] = mjulian2greg(t(i));

    % Cell array for the files to download.
    files = cell(0);

    % Do 2010 first because it straddles the two directories.
    if year == 2010
        if month < 11 && day < 4
            % Use the am data
            prefix = 'am';
            URL = sprintf('%sna/%s/%04d/%02d/%02d', basePath, ...
                prefix, ...
                year, ...
                month, ...
                day);
            for f = 1:ns
                files{f} = sprintf('na%s%04d%02d%02d%02d_%05d_00.pp', ...
                    prefix, ...
                    year, ...
                    month, ...
                    day, ...
                    hour, ...
                    stash(f));
            end
        elseif month > 11 && day > 3
            % Use the mn data
            prefix = 'mn';
            URL = sprintf('%sna/%s/%04d/%02d/%02d', basePath, ...
                prefix, ...
                year, ...
                month, ...
                day);
            files = sprintf('%s_%04d%02d%02d%02d_s00.pp', ...
                prefix, ...
                year, ...
                month, ...
                day, ...
                hour);
        end

    % Check the 2006 data are from the 7th November onwards.
    elseif year == 2006
        if month < 11 
            if day < 7
                error('The MetUM repository does not contain data earlier than 7th November, 2006')
            else
                prefix = 'am';
                URL = sprintf('%sna/%s/%04d/%02d/%02d', basePath, ...
                    prefix, ...
                    year, ...
                    month, ...
                    day);
                for f = 1:ns
                    files{f} = sprintf('na%s%04d%02d%02d%02d_%05d_00.pp', ...
                        prefix, ...
                        year, ...
                        month, ...
                        day, ...
                        hour, ...
                        stash(f));
                end
            end
        end

    % Check the 2012 data are from before the 17th January, 2012.
    elseif year == 2012
        if month > 1
            error('The MetUM repository does not contain data later than 17th January, 2012')
        elseif month == 1
            if day > 17
                error('The MetUM repository does not contain data later than 17th January, 2012')
            else
                prefix = 'mn';
            URL = sprintf('%sna/%s/%04d/%02d/%02d', basePath, ...
                    prefix, ...
                    year, ...
                    month, ...
                    day);
                files = sprintf('%s_%04d%02d%02d%02d_s00.pp', ...
                    prefix, ...
                    year, ...
                    month, ...
                    day, ...
                    hour);
            end
        end

    % Pre-2010 files.
    elseif year < 2010
        % Use the am data.
        prefix = 'am';
            URL = sprintf('%sna/%s/%04d/%02d/%02d', basePath, ...
            prefix, ...
            year, ...
            month, ...
            day);
            for f = 1:ns
                files{f} = sprintf('na%s%04d%02d%02d%02d_%05d_00.pp', ...
                    prefix, ...
                    year, ...
                    month, ...
                    day, ...
                    hour, ...
                    stash(f));
            end

    % Post-2010 files.
    elseif year > 2010
        % Use the mn data.
        prefix = 'mn';
            URL = sprintf('%sna/%s/%04d/%02d/%02d', basePath, ...
            prefix, ...
            year, ...
            month, ...
            day);
        files = sprintf('%s_%04d%02d%02d%02d_s00.pp', ...
            prefix, ...
            year, ...
            month, ...
            day, ...
            hour);
    end

    fprintf('%s: %s\n', URL, files{1})
end

% Close the connection to the FTP server.
close(remote)

if ftbverbose
    fprintf('end   : %s \n', subname)
end


function ftpdata = get_badc_data(remote, URL, files)
% Child function to do the actual downloading from the BADC site via FTP.
% 
% Inputs:
% 
%   remote - FTP object
%   URL - path to the files to download
%   files - cell array of a file or files to download
% 
% Outputs:
% 
%   noidea...

if ~iscell(files)
    error('Provide a cell array of files to download')
end

cd(remote, URL);
nf = length(files);
for i = 1:nf
    tmpdata = mget(remote, files{i});
    ftpdata.(vars{i}).data = 
end


function pp2nc(file, convsh)
% Child function to call the convsh program to convert the obscure pp
% format to a sensible NetCDF which we can more easily read.

% Assume convsh is in /usr/local unless otherwise told.
if nargin == 1
    convsh = '/usr/local/bin/convsh';
end

if exist(file, 'file') ~= 2
    error('File %s not found', file)
end

[path, name, ext] = fileparts(file);
out = fullfile(path, [name, '.nc']);

system([convsh, '-i ', 
