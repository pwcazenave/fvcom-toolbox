function era = read_ERA_wind(Mobj, startDate, endDate, datadir, varlist)
% Reads in ERA Interim files and outputs a struct containing the requested
% variables for a given spatial domain and time period. 
% 
% ERA = read_ERA_wind(Mobj, startDate, endDate, datadir, varlist)
% 
% DESCRIPTION:
%   Find all the ERA Interim NetCDF files in datadir and aggregate them
%   into a single MATLAB struct containing the variables specified in
%   varlist. In addition to the specified variables, time, longitude and
%   latitude will be extracted.
% 
%   The ERA data consist of two files types, gafs*.nc and ggas*.nc. The
%   former contain:
%       - SSTK: sea surface temperature
%       - SP:   surface pressure
%       - MSL:  mean sea level pressure
%       - U10:  u wind (10m)
%       - V10:  v wind (10m)
%       - T2:   2m temperature
%       - SKT:  skin temperature
%   the latter:
%       - UVB:  downward UV radiation at the surface
%       - E:    evaporation
%       - FISR: top of atmosphere incident solar radiation
%       - LSP:  stratiform precipitation (large scale precipitation)
%       - CP:   convective preciptation
%       - SSHF: surface sensible heat flux
%       - SLHF: surface latent heat flux
%       - SSRD: surface solar radiation downwards
%       - SSR:  surface solar radiation
%       - TSR:  top solar radiation
%       - TTR:  top thermal radiation
%       - TP:   total precipitation
%       - TSRC: top net solar radiation, clear sky
%       - TTRC: top upward thermal radiation, clear sky
%       - SSRC: surface net solar radiation, clear sky
%       - STRC: surface net thermal radiation, clear sky
% 
%   Not all of these are necessarily useful with FVCOM, but they're the
%   ones that stood out as potentially useful.
% 
% INPUT:
%   Mobj - Mesh object containing the following fields:
%       * nativeCoords - 'spherical' or 'cartesian'
%       * lon, lat or x, y - node coordinates (depending on nativeCoords)
%   startDate - model start date ([YYYY, MM, DD, hh, mm, ss])
%   endDate -  model end date ([YYYY, MM, DD, hh, mm, ss]). Start and end
%       dates must fall within the same year.
%   datadir - path to the directory which contains the ERA NetCDF files
%   varlist - list of the particular variables to extract from the NetCDF
%       files.
% 
% OUTPUT:
%   era - struct containing the variables specified in varlist.
% 
% Author(s)
%   Pierre Cazenave (Plymouth Marine Laboratory)
% 
% Revision history:
%   2012-10-19 First version based loosely on read_NCEP_wind.m in the
%   fvcom-toolbox.
%   2013-07-08 Modified to automatically find the necessary files and
%   extract the data to cover the time period in question.
% 
%==========================================================================

datadir = '/data/modellers/to_archive/momm-ERA40-interim/';
% datadir = '/users/modellers/pica/Data/ECMWF/2006';
varlist = {'U10', 'V10', 'SKT', 'E', 'TP', 'SSRC', 'STRC'};

%warning off

if nargin ~= 5
    error('Incorrect number of input arguments')
end

subname = 'read_ERA_wind';

global ftbverbose
if ftbverbose
	fprintf('\nbegin : %s \n', subname)
    c = 0; % set a counter for the file load loop
end

if exist(datadir, 'dir') ~= 7
	error('file: %s does not exist', datadir)
end

%--------------------------------------------------------------------------
% Open ERA Interim data and check for time range
%--------------------------------------------------------------------------

[start_year, start_month, start_day] = deal(startDate(1), startDate(2), startDate(3));
[end_year, end_month, end_day] = deal(endDate(1), endDate(2), endDate(3));

if start_year ~= end_year
    error('Cannot (yet?) process more than a single year at a time')
end

% Get a list of files to use. This is not pretty, but seems to work OK.
files = cell(0);
months = start_month:end_month;
for mm = 1:length(months)
    tfiles = dir(fullfile(datadir, num2str(start_year), sprintf('%02d', months(mm)), '*.nc'));
    
    % Add the files to the overall list.
    for f = 1:length(tfiles)
        if strcmpi(tfiles(f).name, '.') || strcmpi(tfiles(f).name, '..')
            continue
        else
            files{length(files) + 1} = fullfile(datadir, num2str(start_year), sprintf('%02d', months(mm)), tfiles(f).name);
        end
    end
end

nf = length(files);

% Get the times
tggas = [];
tgafs = [];
gg = 0;
ga = 0;
for f = 1:nf
    
    % Get the filename only for prefix comparison
    [~, ff, ext] = fileparts(files{f});

    if ftbverbose
        c = c + 1;
        fprintf('File %s (%i of %i)... ', [ff, '.', ext], c, nf)
    end

    % File name is ????YYYYMMDDOOFF or ????YYYYMMDDHHFF (with no apparent
    % rhyme or reason...). OO is the time origin (midday/midnight) and FF
    % is the number of hours since the time origin (3, 6, 9 or 12). Why the
    % files couldn't consistently use 24 hour times...
    ymd = greg2mjulian(str2double(ff(5:8)), ...
        str2double(ff(9:10)), ...
        str2double(ff(11:12)), ...
        str2double(ff(13:14)) + str2double(ff(15:16)), ...
        0, 0); % don't have minutes and seconds in the file name

    if strcmpi(ff(1:4), 'ggas')
        gg = gg + 1;
        tggas(gg) = ymd;
    elseif strcmpi(ff(1:4), 'gafs')
        ga = ga + 1;
        tgafs(ga) = ymd;
    else
        warning('Unrecognised ERA Interim file prefix (%s)', ff(1:4))
    end
    
    % Now load in the variables in varlist for all the relevant files into two
    % structs, ggas and gafs. We're assuming that the files are listed
    % increasing with time (i.e. ????YYYYMMDDhhmm.nc).
    nc = netcdf.open(files{f}, 'NOWRITE');

    % Only do the spatial data on the first file (the data are global).
    if f == 1
        lat_varid = netcdf.inqVarID(nc, 'latitude');
        lon_varid = netcdf.inqVarID(nc, 'longitude');
        eralatvector = netcdf.getVar(nc, lat_varid);
        eralonvector = netcdf.getVar(nc, lon_varid);
        [era.lon, era.lat] = meshgrid(eralonvector, eralatvector);

        % To save on memory, cull the data which doesn't come from the
        % domain we've specified (with a 2 element buffer).
        if min(Mobj.lon) < 0 && max(Mobj.lon) < 0
            % Easy, just shift by 360.
            idx_lon = find(eralonvector > min(Mobj.lon) + 360 & eralonvector < max(Mobj.lon) + 360);
        elseif min(Mobj.lon) < 0 && max(Mobj.lon) > 0
            % This is the tricky one. We'll do two passes to extract
            % the western chunk first (min(Mobj.lon) + 360 to 360),
            % then the eastern chunk (0 - max(Mobj.lon))
            idx_lon{1} = find(eralonvector >= min(Mobj.lon) + 360);
            idx_lon{2} = find(eralonvector < max(Mobj.lon));
        else
            % Dead easy, we're in the eastern hemisphere, so nothing
            % too strenuous here.
            idx_lon = find(eralonvector > min(Mobj.lon) & eralonvector < max(Mobj.lon));
        end
        % Latitudes are easy because there's only one system.
        idx_lat = find(eralatvector > min(Mobj.lat) & eralatvector < max(Mobj.lat));
    end

    if f == 1
        era.time = tstart;
    else
        era.time = [era.time; tstart];
    end

    for v = 1:length(varlist)
        % Use a try catch to pass on the files which don't contain the
        % variable we're currently looping over.
        try
            % Get the data
            varid_era = netcdf.inqVarID(nc, varlist{v});
            data = netcdf.getVar(nc, varid_era, 'single');
        catch
            % Skip this variable and try the next one
            continue
        end

        try
            % Try to apply the scale factor and offset.
            sf = netcdf.getAtt(nc, varid_era, 'scale_factor', 'double');
            ao = netcdf.getAtt(nc, varid_era, 'add_offset', 'double');
            if exist(sprintf('era.(''%s'').data', varlist{v}), 'var') == 1
                era.(varlist{v}).data = permute(double(ao + (data(idx_lon, idx_lat) .* sf)), [2, 1, 3]);
            else
                era.(varlist{v}).data = cat(3, era.(varlist{v}).data, permute(double(ao + (data(idx_lon, idx_lat) .* sf)), [2, 1, 3]));
            end
        catch
            % Otherwise just dump the data as is.
            if exist(sprintf('era.(''%s'').data', varlist{v}), 'var') == 1
                era.(varlist{v}).data = permute(double(data(idx_lon, idx_lat)), [2, 1, 3]);
            else
                era.(varlist{v}).data = cat(3, era.(varlist{v}).data, permute(double(data(idx_lon, idx_lat)), [2, 1, 3]));
            end
        end

        % Put the lon/lat into the field too.
        era.(varlist{v}).lat = eralatvector(idx_lat);
        era.(varlist{v}).lon = eralonvector(idx_lon);

    end
    
    if ftbverbose
        fprintf('done.\n')
    end
end

netcdf.close(nc)

if ftbverbose
    fprintf('end   : %s \n', subname)
end
