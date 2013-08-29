function met_files = get_MetUM_pp(modelTime, credentials)
% Get the required parameters from the Met Office Unified Model (TM)
% (hereafter MetUM) to use in FVCOM surface forcing.
%
% met_files = get_MetUM_pp(modelTime, credentials)
%
% DESCRIPTION:
%   Using FTP access, extract the necessary parameters to create an FVCOM
%   forcing file. Requires the air_sea toolbox (see below for where to get
%   it). Data are sampled four times daily.
%
% INPUT:
%   modelTime - Modified Julian Date start and end times array
%   credentials - struct with fields username and password to access the
%   FTP server.
%
% OUTPUT:
%   met_files - cell array of file names downloaded from the BADC servers.
%
% The PP files downloaded give:
%     - surface_net_downward_shortwave_flux (W m-2)
%     - surface_downwelling_shortwave_flux_in_air (W m-2)
%     - surface_net_downward_longwave_flux (W m-2)
%     - surface_downwelling_longwave_flux_in_air (W m-2)
%     - surface_upward_latent_heat_flux (W m-2)
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
% EXAMPLE USAGE:
%   met_files = get_MetUM_pp([51725, 51757], {'username', 'password'});
%
% TODO:
%   Add support for the AP directories on the FTP server.
%
% Author(s)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2013-05-07 First version.
%   2013-06-24 Update some of the code to correctly separate the different
%   directories (e.g. am vs mn in 2010). Also farm out the FTP and
%   conversion from PP format to NetCDF to separate functions
%   (get_BADC_data.m and pp2nc.m, respectively). Renamed the function from
%   get_MetUM_forcing to get_MetUM_pp to better reflect what it does.
%
%==========================================================================

subname = 'get_MetUM_forcing';

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

nt = ceil(modelTime(end)) - floor(modelTime(1));
if nt > 365
    error('Can''t (yet) process more than a year at a time.')
end

yearStart = mjulian2greg(modelTime(1));
yearEnd = mjulian2greg(modelTime(end));

% Four times daily outputs at 0000, 0600, 1200 and 1800
t = modelTime(1):1/4:modelTime(end);

assert(yearEnd == yearStart, 'Can''t (yet) process across a year boundary.')
assert(yearStart >= 2006 && yearEnd <= 2012, 'The MetUM repository does not contain data earlier than 2006 and later than 2012')

% For the pre-2010 data, we need to download several files with
% unique names. The names are based on the STASH numbers and the date:
%   naamYYYYMMDDHH_STASH#_00.pp
% The numbers we're interested in are stored in stash.
stash = [2, 3, 407, 408, 409, 4222, 9229, 16004, ...
    1201, 1235, 2207, 2201, 3217, 3225, 3226, 3234, 3236, 3237, 3245, ...
    5216, 16222, 20004];
% The stash numbers and their corresponding forcing type.
%
% AP = analysis, pressure levels
% AM = analysis, model levels
%
% |---------|-------------------------------------------|
% | stash # | forcing type                              |
% |-----------------------------------------------------|
% |                         AM                          |
% |-----------------------------------------------------|
% | 2       | eastward_wind / x_wind                    |
% | 3       | northward_wind / y_wind                   |
% | 407     | air_pressure                              |
% | 408     | air_pressure                              |
% | 409     | surface_air_pressure                      |
% | 4222    | RAILFALL RATE OUT OF MODEL LEVELS         |
% | 9229    | RELATIVE HUMIDITY AFTER MAIN CLOUD        |
% | 16004   | air_temperature                           |
% |-----------------------------------------------------|
% |                         AP                          |
% |-----------------------------------------------------|
% | 1201    | surface_net_downward_shortwave_flux       |
% | 1235    | surface_downwelling_shortwave_flux_in_air |
% | 2207    | surface_downwelling_longwave_flux_in_air  |
% | 2201    | surface_net_downward_longwave_flux        |
% | 3217    | surface_upward_sensible_heat_flux         |
% | 3225    | eastward_wind / x_wind                    |
% | 3226    | northward_wind / y_wind                   |
% | 3234    | surface_upward_latent_heat_flux           |
% | 3236    | air_temperature                           |
% | 3237    | specific_humidity                         |
% | 3245    | relative_humidity                         |
% | 5216    | precipitation_flux                        |
% | 16222   | air_pressure_at_sea_level                 |
% | 20004   | [RIVER OUTFLOW]                           |
% |---------|-------------------------------------------|
%

ns = length(stash);

% From where will we be downloading the data?
site = 'ftp.ceda.ac.uk';
basePath = 'badc/ukmo-um/data/nae/';

% Depending on the year we're extracting, we need to append different
% directories to get the data. 
for i = 1:nt * 4 % four files per day (at 0000, 0600, 1200 and 1800).

    [year, month, day, hour] = mjulian2greg(t(i));

    % Cell array for the files to download.
    files = cell(0);

    % Do 2010 first because it straddles the two directories.
    if year == 2010
        % Use modified julian dates for the thresholds for each directory.
        amthresh = greg2mjulian(2010, 11, 03, 00, 00, 00);

        if greg2mjulian(year, month, day, hour, 00, 00) <= amthresh
            % Use the am data
            prefix = 'am';
            filepath = sprintf('%sna/%s/%04d/%02d/%02d', basePath, ...
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
        else
            % Use the mn data
            prefix = 'mn';
            filepath = sprintf('%sna/%s/%04d/%02d/%02d', basePath, ...
                prefix, ...
                year, ...
                month, ...
                day);
            files = {sprintf('%s_%04d%02d%02d%02d_s00.pp', ...
                prefix, ...
                year, ...
                month, ...
                day, ...
                hour)};
        end
        sprintf('%s', filepath);

    % Check the 2006 data are from the 7th November onwards.
    elseif year == 2006
        if month < 11 
            if day < 7
                error('The MetUM repository does not contain data earlier than 7th November, 2006')
            else
                prefix = 'am';
                filepath = sprintf('%sna/%s/%04d/%02d/%02d', basePath, ...
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
                filepath = sprintf('%sna/%s/%04d/%02d/%02d', basePath, ...
                    prefix, ...
                    year, ...
                    month, ...
                    day);
                files = {sprintf('%s_%04d%02d%02d%02d_s00.pp', ...
                    prefix, ...
                    year, ...
                    month, ...
                    day, ...
                    hour)};
            end
        end

    % Pre-2010 files.
    elseif year < 2010
        % Use the am data.
        prefix = 'am';
            filepath = sprintf('%sna/%s/%04d/%02d/%02d', basePath, ...
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
        % Use the sn data (has everything we need but doesn't have 70
        % vertical levels!).
        prefix = 'sn';
        filepath = sprintf('%sna/%s/%04d/%02d/%02d', basePath, ...
            prefix, ...
            year, ...
            month, ...
            day);
        files = {sprintf('%s_%04d%02d%02d%02d_s00.pp', ...
            prefix, ...
            year, ...
            month, ...
            day, ...
            hour)};
    end

    met_files{i} = get_BADC_data(site, filepath, files, credentials);

end

if ftbverbose
    fprintf('end   : %s \n', subname)
end
