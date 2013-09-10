function Mobj = interp_POLCOMS2FVCOM(Mobj, ts, start_date, varlist)
% Use an FVCOM restart file to seed a model run with spatially varying
% versions of otherwise constant variables.
%
% function interp_POLCOMS2FVCOM(Mobj, ts, start_date, varlist)
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
%                   - Mobj.lon, Mobj.lat - node coordinates (long/lat)
%                   - Mobj.lonc, Mobj.latc - element coordinates (long/lat)
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
%   interp_POLCOMS2FVCOM(Mobj, '/tmp/ts.nc', [2006, 01, 01, 00, 00, 00], ...
%       {'lon', 'lat', 'ETWD', 'x1XD', 'ucurD', 'vcurD', 'rholocalD', 'time'})
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2013-02-08 First version.
%   2013-05-16 Add support for parallel for-loops (not mandatory, but
%   enabled if the Parallel Computing Toolbox is available).
%   2013-06-06 Fix the vertical ordering of the POLCOMS data. POLCOMS'
%   scalar values (temperature, salinity etc.) are stored seabed to
%   surface; its depths are stored surface to seabed; FVCOM stores
%   everything surface to seabed. As such, the POLCOMS scalar values need
%   to be flipped upside down to match everything else.
%   2013-07-30 Add density and u/v velocity components to the variables to
%   interpolate.
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
fe = numel(Mobj.lonc);

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
density = permute(squeeze(pc.rholocalD.data(:, :, :, tidx)), [2, 1, 3]);
u = permute(squeeze(pc.ucurD.data(:, :, :, tidx)), [2, 1, 3]);
v = permute(squeeze(pc.vcurD.data(:, :, :, tidx)), [2, 1, 3]);
depth = permute(squeeze(pc.depth.data(:, :, :, tidx)), [2, 1, 3]);
mask = depth(:, :, end) >= 0; % land is positive.

pc.tempz = grid_vert_interp(Mobj, lon, lat, temperature, depth, mask);
pc.salz = grid_vert_interp(Mobj, lon, lat, salinity, depth, mask);
pc.denz = grid_vert_interp(Mobj, lon, lat, density, depth, mask);
pc.uvelz = grid_vert_interp(Mobj, lon, lat, u, depth, mask);
pc.vvelz = grid_vert_interp(Mobj, lon, lat, v, depth, mask);

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
fvdens = nan(fn, fz);
fvuvel = nan(fe, fz);
fvvvel = nan(fe, fz);

plon = lon(:);
plat = lat(:);
flon = Mobj.lon;
flat = Mobj.lat;
flonc = Mobj.lonc;
flatc = Mobj.latc;
ptempz = pc.tempz;
psalz = pc.salz;
pdenz = pc.denz;
puvelz = pc.uvelz;
pvvelz = pc.vvelz;

tic
parfor zi = 1:fz
    % Set up the interpolation objects.
    ft = TriScatteredInterp(plon, plat, reshape(ptempz(:, :, zi), [], 1), 'natural');
    fs = TriScatteredInterp(plon, plat, reshape(psalz(:, :, zi), [], 1), 'natural');
    fd = TriScatteredInterp(plon, plat, reshape(pdenz(:, :, zi), [], 1), 'natural');
    fu = TriScatteredInterp(plon, plat, reshape(puvelz(:, :, zi), [], 1), 'natural');
    fv = TriScatteredInterp(plon, plat, reshape(pvvelz(:, :, zi), [], 1), 'natural');
    % Interpolate temperature and salinity onto the unstructured grid.
    fvtemp(:, zi) = ft(flon, flat);
    fvsalt(:, zi) = fs(flon, flat);
    fvdens(:, zi) = fd(flon, flat);
    fvuvel(:, zi) = fu(flonc, flatc);
    fvvvel(:, zi) = fv(flonc, flatc);
end

clear plon plat flon flat flonc flatc ptempz psalz pdenz puvelz pvvelz

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
    fvdens(fvnanidx(ni), :) = fvdens(fvfinidx(di), :);
    fvuvel(fvnanidx(ni), :) = fvuvel(fvfinidx(di), :);
    fvvvel(fvnanidx(ni), :) = fvvvel(fvfinidx(di), :);
end

if ftbverbose
    fprintf('done.\n') 
    toc
end

Mobj.restart.temp = fvtemp;
Mobj.restart.salinity = fvsalt;
Mobj.restart.rho1 = fvdens;
Mobj.restart.u = fvuvel;
Mobj.restart.v = fvvvel;

% Close the MATLAB pool if we opened it.
if wasOpened
    matlabpool close
end

if ftbverbose
    fprintf('end   : %s\n', subname)
end

%% Debugging figure

% close all
%
% ri = 85; % column index
% ci = 95; % row index
%
% [~, idx] = min(sqrt((Mobj.lon - lon(ri, ci)).^2 + (Mobj.lat - lat(ri, ci)).^2));
%
% % Vertical profiles
% figure
% clf
%
% % The top row shows the temperature/salinity values as plotted against
% % index (i.e. position in the array). I had thought POLCOMS stored its data
% % from seabed to sea surface, but looking at the NetCDF files in ncview,
% % it's clear that the data are in fact stored surface to seabed (like
% % FVCOM). As such, the two plots in the top row should be upside down (i.e.
% % surface values at the bottom of the plot). The two lower rows should have
% % three lines which all match: the raw POLCOMS data, the POLCOMS data for
% % the current time step (i.e. that in 'temperature' and 'salinity') and the
% % interpolated FVCOM data against depth.
% %
% % Also worth noting, the pc.*.data have the rows and columns flipped, so
% % (ci, ri) in pc.*.data and (ri, ci) in 'temperature', 'salinity' and
% % 'depth'. Needless to say, the two lines in the lower plots should
% % overlap.
%
% subplot(2,2,1)
% plot(squeeze(pc.ETWD.data(ci, ri, :, tidx)), 1:size(depth, 3), 'rx:')
% xlabel('Temperature (^{\circ}C)')
% ylabel('Array index')
% title('Array Temperature')
%
% subplot(2,2,2)
% plot(squeeze(pc.x1XD.data(ci, ri, :, tidx)), 1:size(depth, 3), 'rx:')
% xlabel('Salinity')
% ylabel('Array index')
% title('Array Salinity')
%
% subplot(2,2,3)
% % Although POLCOMS stores its temperature values from seabed to surface,
% % the depths are stored surface to seabed. Nice. Flip the
% % temperature/salinity data accordingly. The lines here should match one
% % another.
% plot(squeeze(pc.ETWD.data(ci, ri, :, tidx)), squeeze(pc.depth.data(ci, ri, :, tidx)), 'rx-')
% hold on
% plot(squeeze(temperature(ri, ci, :)), squeeze(depth(ri, ci, :)), '.:')
% % Add the equivalent interpolated FVCOM data point
% plot(fvtemp(idx, :), Mobj.siglayz(idx, :), 'g.-')
% xlabel('Temperature (^{\circ}C)')
% ylabel('Depth (m)')
% title('Depth Temperature')
% legend('pc', 'temp', 'fvcom', 'location', 'north')
% legend('boxoff')
%
% subplot(2,2,4)
% % Although POLCOMS stores its temperature values from seabed to surface,
% % the depths are stored surface to seabed. Nice. Flip the
% % temperature/salinity data accordingly. The lines here should match one
% % another.
% plot(squeeze(pc.x1XD.data(ci, ri, :, tidx)), squeeze(pc.depth.data(ci, ri, :, tidx)), 'rx-')
% hold on
% plot(squeeze(salinity(ri, ci, :)), squeeze(depth(ri, ci, :)), '.:')
% % Add the equivalent interpolated FVCOM data point
% plot(fvsalt(idx, :), Mobj.siglayz(idx, :), 'g.-')
% xlabel('Salinity')
% ylabel('Depth (m)')
% title('Depth Salinity')
% legend('pc', 'salt', 'fvcom', 'location', 'north')
% legend('boxoff')
%
% %% Plot the sample location
% figure
% dx = mean(diff(pc.lon.data));
% dy = mean(diff(pc.lat.data));
% z = depth(:, :, end); % water depth (bottom layer depth)
% z(mask) = 0; % clear out nonsense values
% pcolor(lon - (dx / 2), lat - (dy / 2), z)
% shading flat
% axis('equal', 'tight')
% daspect([1.5, 1, 1])
% colorbar
% caxis([-150, 0])
% hold on
% plot(lon(ri, ci), lat(ri, ci), 'ko', 'MarkerFaceColor', 'w')

