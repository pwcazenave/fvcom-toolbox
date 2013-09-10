function Mobj = interp_HYCOM2FVCOM(Mobj, hycom, start_date, varlist)
% Use an FVCOM restart file to seed a model run with spatially varying
% versions of otherwise constant variables.
%
% function interp_HYCOM2FVCOM(Mobj, hycom, start_date, varlist)
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
%   hycom       = Struct output by get_HYCOM_forcing. Must include fields:
%                   - hycom.lon, hycom.lat - rectangular arrays.
%                   - hycom.Depth - HYCOM depth levels.
%   start_date  = Gregorian start date array (YYYY, MM, DD, hh, mm, ss).
%   varlist     = cell array of fields to use from the HYCOM struct.
%
% OUTPUT:
%   Mobj.restart = struct whose field names are the variables which have
%   been interpolated (e.g. Mobj.restart.temperature for HYCOM
%   temperature).
%
% EXAMPLE USAGE
%   interp_HYCOM2FVCOM(Mobj, hycom, [2006, 01, 01, 00, 00, 00], ...
%       {'lon', 'lat', 'time', 'temperature', 'salinity'})
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2013-09-10 First version based on interp_POLCOMS2FVCOM.m.
%
%==========================================================================

subname = 'interp_HYCOM2FVCOM';

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

% Number of sigma layers, grid nodes and elements.
[fn, fz] = size(Mobj.siglayz);
fe = numel(Mobj.lonc);

% Given our intput time (in start_date), find the nearest time index for
% the HYCOM data.
stime = greg2mjulian(start_date(1), start_date(2), ...
    start_date(3), start_date(4), ...
    start_date(5), start_date(6));
[~, tidx] = min(abs(Mobj.ts_times - stime));

for vv = 1:length(varlist);
    
    currvar = varlist{vv};

    switch currvar
        case {'lon', 'lat', 'longitude', 'latitude', 't_1', 'time', 'Depth'}
            continue

        otherwise

            %--------------------------------------------------------------
            % Interpolate the regularly gridded data onto the FVCOM grid
            % (vertical grid first).
            %--------------------------------------------------------------

            if ftbverbose
                fprintf('%s : interpolate %s vertically... ', subname, currvar)
            end

            [hx, hy, ~, ~] = size(hycom.(currvar).data);
            hdepth = permute(repmat(-hycom.Depth.data, [1, hx, hy]), [2, 3, 1]);

            % Make a land mask of the HYCOM domain (based on the surface
            % layer from the first time step).
            landmask = hycom.(currvar).data(:, :, 1, 1) > 1.26e29;

            % Create a mask of all the layers which are deeper than the
            % water depth.
            zmask = hycom.(currvar).data(:, :, :, 1) > 1.26e29;

            % Set invalid depths to NaN.
            hdepth(zmask) = nan;

            ftbverbose = false;
            hyinterp.(currvar) = grid_vert_interp(Mobj, ...
                hycom.lon, hycom.lat, ...
                squeeze(hycom.(currvar).data(:, :, :, tidx)), ...
                hdepth, landmask);
            ftbverbose = true;

            % Clear out additional values which have fallen through the
            % mask check above.
            hyinterp.(currvar)(hyinterp.(currvar) > 1.26e20) = nan;

            %--------------------------------------------------------------
            % Now we have vertically interpolated data, we can interpolate
            % each sigma layer onto the FVCOM unstructured grid ready to
            % write out to NetCDF. We'll use the triangular interpolation
            % in MATLAB with the natural method (gives pretty good results,
            % at least qualitatively).
            %--------------------------------------------------------------

            if ftbverbose
                fprintf('horizontally... ')
            end

            fvtemp = nan(fn, fz);
            if strcmpi(currvar, 'u') || strcmpi(currvar, 'v')
                fvtemp = nan(fe, fz);
            end
            plon = hycom.lon(:);
            plat = hycom.lat(:);
            flon = Mobj.lon;
            flat = Mobj.lat;
            flonc = Mobj.lonc;
            flatc = Mobj.latc;
            hytempz = hyinterp.(currvar);

            tic
            parfor zi = 1:fz
                % Set up the interpolation object and interpolate the
                % current variable to the FVCOM unstructured grid.
                ft = TriScatteredInterp(plon, plat, reshape(hytempz(:, :, zi), [], 1), 'natural');
                if strcmpi(currvar, 'u') || strcmpi(currvar, 'v')
                    fvtemp(:, zi) = ft(flonc, flatc);
                else
                    fvtemp(:, zi) = ft(flon, flat);
                end
            end

            % Unfortunately, TriScatteredInterp won't extrapolate,
            % returning instead NaNs outside the original data's extents.
            % So, for each NaN position, find the nearest non-NaN value and
            % use that instead. The order in which the NaN-nodes are found
            % will determine the spatial pattern of the extrapolation.

            % We can assume that all layers will have NaNs in the same
            % place (horizontally), so just use the surface layer (1) for
            % the identification of NaNs. Also store the finite values so
            % we can find the nearest real value to the current NaN node
            % and use its temperature and salinity values.
            fvidx = 1:fn;
            fvnanidx = fvidx(isnan(fvtemp(:, 1)));
            fvfinidx = fvidx(~isnan(fvtemp(:, 1)));

            for ni = 1:length(fvnanidx)
                xx = Mobj.lon(fvnanidx(ni));
                yy = Mobj.lat(fvnanidx(ni));
                % Find the nearest non-nan temperature and salinity value.
                [~, di] = min(sqrt((Mobj.lon(fvfinidx) - xx).^2 + (Mobj.lat(fvfinidx) - yy).^2));

                fvtemp(fvnanidx(ni), :) = fvtemp(fvfinidx(di), :);

            end

            clear plon plat flon flat flonc flatc ptempz

            Mobj.restart.(currvar) = fvtemp;

            te = toc;

            if ftbverbose
                fprintf('done. (elapsed time = %.2f seconds)\n', te) 
            end

    end
end

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

