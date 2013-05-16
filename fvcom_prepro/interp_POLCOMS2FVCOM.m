function Mobj = interp_POLCOMS2FVCOM(Mobj, ts, start_date, varlist)
% Use an FVCOM restart file to seed a model run with spatially varying
% versions of otherwise constant variables (temperature and salinity only
% for the time being).
%
% function interp_POLCOMS2FVCOM(Mobj, ts, start_date, fv_restart, varlist)
%
% DESCRIPTION:
%    FVCOM does not yet support spatially varying temperature and salinity
%    inputs as initial conditions. To avoid having to run a model for a
%    long time in order for temperature and salinity to settle within the
%    model from the atmospheric and boundary forcing, we can use a restart
%    file to cheat. For this, we need temperature and salinity
%    (potentially other variables too) interpolated onto the unstructured
%    grid. The interpolated data can then be written out with
%    write_FVCOM_restart.m.
%
% INPUT:
%   Mobj        = MATLAB mesh structure which must contain:
%                   - Mobj.siglayz - sigma layer depths for all model
%                   nodes.
%                   - Mobj.lon, Mobj.lat - node coordinates (long/lat).
%   ts          = Cell array of POLCOMS AMM NetCDF file(s) in which 4D
%   variables of temperature and salinity (called 'ETWD' and 'x1XD') exist.
%   Its/their shape should be (y, x, sigma, time).
%   start_date  = Gregorian start date array (YYYY, MM, DD, hh, mm, ss).
%   varlist     = cell array of variables to extract from the NetCDF files.
% 
% OUTPUT:
%   Mobj.restart = struct whose field names are the variables which have
%   been interpolated (e.g. Mobj.restart.ETWD for POLCOMS daily mean
%   temperature).
%
% EXAMPLE USAGE
%   interp_POLCOMS2FVCOM(Mobj, '/tmp/ts.nc', '2006-01-01 00:00:00', ...
%       {'lon', 'lat', 'ETWD', 'x1XD', 'time'})
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2013-02-08 First version.
%   2013-05-16 Add support for parallel for-loops (not mandatory, but
%   enabled if the Parallel Computing Toolbox is available).
%
%==========================================================================

subname = 'interp_POLCOMS2FVCOM';

global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

% Run jobs on multiple workers if we have that functionality. Not sure if
% it's necessary, but check we have the Parallel Toolbox first.
wasOpened = false;
if license('test', 'Distrib_Computing_Toolbox')
    % We have the Parallel Computing Toolbox, so launch a bunch of workers.
    if matlabpool('size') == 0
        % Force pool to be local in case we have remote pools available.
        matlabpool open local
        wasOpened = true;
    end
end

%--------------------------------------------------------------------------
% Extract the NetCDF data specified in varlist
%--------------------------------------------------------------------------

% Data format:
% 
%   pc.ETWD.data and pc.x1XD.data are y, x, sigma, time
% 
pc = get_POLCOMS_netCDF(ts, varlist);

% Number of sigma layers.
[fn, fz] = size(Mobj.siglayz);

% Make rectangular arrays for the nearest point lookup.
[lon, lat] = meshgrid(pc.lon.data, pc.lat.data);

% Convert the current times to Modified Julian Day (this is a bit ugly).
pc.time.all = strtrim(regexp(pc.time.units, 'since', 'split'));
pc.time.datetime = strtrim(regexp(pc.time.all{end}, ' ', 'split'));
pc.time.ymd = str2double(strtrim(regexp(pc.time.datetime{1}, '-', 'split')));
pc.time.hms = str2double(strtrim(regexp(pc.time.datetime{2}, ':', 'split')));

Mobj.ts_times = greg2mjulian(...
    pc.time.ymd(1), ...
    pc.time.ymd(2), ...
    pc.time.ymd(3), ...
    pc.time.hms(1), ...
    pc.time.hms(2), ...
    pc.time.hms(3)) + (pc.time.data / 3600 / 24);

% Given our intput time (in start_date), find the nearest time
% index for the regularly gridded data.
stime = greg2mjulian(start_date(1), start_date(2), ...
    start_date(3), start_date(4), ...
    start_date(5), start_date(6));
[~, tidx] = min(abs(Mobj.ts_times - stime));

%--------------------------------------------------------------------------
% Interpolate the regularly gridded data onto the FVCOM grid (vertical grid
% first).
%--------------------------------------------------------------------------

if ftbverbose
    fprintf('%s : interpolate POLCOMS onto FVCOM''s vertical grid... ', subname)
end

% Permute the arrays to be x by y rather than y by x.
temperature = permute(squeeze(pc.ETWD.data(:, :, :, tidx)), [2, 1, 3]);
salinity = permute(squeeze(pc.x1XD.data(:, :, :, tidx)), [2, 1, 3]);
depth = permute(squeeze(pc.depth.data(:, :, :, tidx)), [2, 1, 3]);
mask = depth(:, :, end) >= 0; % land is positive.

pc.tempz = grid_vert_interp(Mobj, lon, lat, temperature, depth, mask);
pc.salz = grid_vert_interp(Mobj, lon, lat, salinity, depth, mask);

if ftbverbose
    fprintf('done.\n') 
end

%--------------------------------------------------------------------------
% Now we have vertically interpolated data, we can interpolate each sigma
% layer onto the FVCOM unstructured grid ready to write out to NetCDF.
% We'll use the triangular interpolation in MATLAB with the natural method
% (gives pretty good results, at least qualitatively).
%--------------------------------------------------------------------------

if ftbverbose
    fprintf('%s : interpolate POLCOMS onto FVCOM''s horizontal grid... ', subname)
end

fvtemp = nan(fn, fz);
fvsalt = nan(fn, fz);

tic
parfor zi = 1:fz
    % Set up the interpolation objects.
    ft = TriScatteredInterp(lon(:), lat(:), reshape(pc.tempz(:, :, zi), [], 1), 'natural');
    fs = TriScatteredInterp(lon(:), lat(:), reshape(pc.salz(:, :, zi), [], 1), 'natural');
    % Interpolate temperature and salinity onto the unstructured grid.
    fvtemp(:, zi) = ft(Mobj.lon, Mobj.lat);
    fvsalt(:, zi) = fs(Mobj.lon, Mobj.lat);
end

% Unfortunately, TriScatteredInterp won't extrapolate, returning instead
% NaNs outside the original data's extents. So, for each NaN position, find
% the nearest non-NaN value and use that instead. The order in which the
% NaN-nodes are found will determine the spatial pattern of the
% extrapolation.

% We can assume that all layers will have NaNs in the same place
% (horizontally), so just use the surface layer (1) for the identification
% of NaNs. Also store the finite values so we can find the nearest real
% value to the current NaN node and use its temperature and salinity
% values.
fvidx = 1:fn;
fvnanidx = fvidx(isnan(fvtemp(:, 1)));
fvfinidx = fvidx(~isnan(fvtemp(:, 1)));

% Can't parallelise this one (easily). It shouldn't be a big part of the
% run time if your source data covers the domain sufficiently.
for ni = 1:length(fvnanidx)
    % Current position
    xx = Mobj.lon(fvnanidx(ni));
    yy = Mobj.lat(fvnanidx(ni));
    % Find the nearest non-nan temperature and salinity value.
    [~, di] = min(sqrt((Mobj.lon(fvfinidx) - xx).^2 + (Mobj.lat(fvfinidx) - yy).^2));
    % Replace the temperature and salinity values at all depths at the
    % current NaN position with the closest non-nan value.
    fvtemp(fvnanidx(ni), :) = fvtemp(fvfinidx(di), :);
    fvsalt(fvnanidx(ni), :) = fvsalt(fvfinidx(di), :);
end

if ftbverbose
    fprintf('done.\n') 
    toc
end

Mobj.restart.temp = fvtemp;
Mobj.restart.salinity = fvsalt;

% Close the MATLAB pool if we opened it.
if wasOpened
    matlabpool close
end

if ftbverbose
    fprintf('end   : %s\n', subname)
end
