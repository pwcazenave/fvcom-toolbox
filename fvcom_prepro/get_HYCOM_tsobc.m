function Mobj = get_HYCOM_tsobc(Mobj, hycom, varlist)
% Read temperature, salinity, u and v data from the HYCOM model output
% structure and interpolate onto the open boundaries in Mobj.
%
% function Mobj = get_HYCOM_tsobc(Mobj, hycom, varlist)
%
% DESCRIPTION:
%    Interpolate temperature and salinity values onto the FVCOM open
%    boundaries at all sigma levels.
%
% INPUT:
%   Mobj    = MATLAB mesh structure which must contain:
%               - Mobj.siglayz - sigma layer depths for all model nodes.
%               - Mobj.siglayzc - sigma layer depths for all model elements.
%               - Mobj.lon, Mobj.lat - node coordinates (lat/long).
%               - Mobj.lonc, Mobj.latc - element coordinates (lat/long).
%               - Mobj.read_obc_nodes - cell array of open boundary nodes.
%               - Mobj.read_obc_elems - cell array of open boundary
%               elements (only if using velocities - use find_nested_region
%               to get these indices).
%               - Mobj.h - water depths at nodes.
%               - Mobj.tri - triangulation table for the grid (nv in FVCOM
%               terms).
%               - Mobj.nObcNodes - number of nodes in each open boundary.
%   hycom   = Struct with HYCOM data covering the model domain. Unless
%             varlist is specified (see below), all 4D fields will be
%             interpolated onto the open boundaries (1-3D data will be
%             ignored).
%   varlist = [optional] cell array of variable (field) names to use from
%             hycom.
%
% OUTPUT:
%    Mobj = MATLAB structure with new fields whose names match those given
%    in hycom. The fields have sizes (sum(Mobj.nObcNodes), sigma, time).
%    The time dimension is determined based on the number of time steps in
%    hycom. The ts_time variable is just the input file times in Modified
%    Julian Day.
%
% EXAMPLE USAGE
%    hycom = get_HYCOM_forcing(Mobj, [51500, 51531]);
%    Mobj = get_HYCOM_tsobc(Mobj, hycom, {'u', 'v', 'temperature', 'salinity'})
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2013-09-03 First version based on get_POLCOMS_tsobc.m.
%    2014-04-28 Update interp1 function to use pchip instead of csap as the
%    latter will be removed in a future version of MATLAB and the
%    innumerable warnings were doing my nut in. I checked the output using
%    the new interp1 call and it's identical to the old version. Also
%    update the parallel toolbox stuff for the same reason (future
%    removal).
%    2015-05-21 Remove the old parallel processing bits and replace with
%    the current versions.
%    2016-03-15 Add fallback interpolation to inverse distance weighted if
%    the triangular interpolation fails (which can happen if the points
%    supplied are all in a line, for example).
%    2017-01-27 Subset the coarse data (HYCOM, CMEMS etc.). This yields a
%    significant performance improvement (both in memory usage and time).
%    2017-02-16 Further performance improvement by only using the coarse
%    data in the vicinity of the open boundary nodes.
%    2017-10-12 Fix bug in indexing the open boundary positions which may
%    have caused the interpolation to fail as the identified positions
%    would be too far from the open boundary nodes.
%
%==========================================================================

[~, subname] = fileparts(mfilename('fullpath'));

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

if license('test', 'Distrib_Computing_Toolbox')
    % We have the Parallel Computing Toolbox, so launch a bunch of workers.
    if isempty(gcp('nocreate'))
        % Force pool to be local in case we have remote pools available.
        parpool('local');
    end
end

if nargin == 2
    fields = fieldnames(hycom);
else
    % We always want Depth because we need it to interpolate the vertical
    % component.
    fields = unique(['Depth', varlist], 'stable');
end

assert(isfield(hycom, 'Depth'), 'Require a depth vector to interpolate vertically.')

% Find the first 4D array and use it to get the number of vertical levels
% and time steps.
for ff = 1:length(fields)
    if isfield(hycom.(fields{ff}), 'data') && ndims(hycom.(fields{ff}).data) > 3
        [nx, ny, nz, nt] = size(hycom.(fields{ff}).data);
        break
    end
end

% Use the existing rectangular arrays for the nearest point lookup.
[lon, lat] = deal(hycom.lon, hycom.lat);

% Find the indices of all the HYCOM data which cover the open boundary
% nodes given by some distance. This should massively increase performance
% by chopping out the vast majority of the data. Find the 6 closest ones
% for each node/element.

obc_hycom_idx = [];
if isfield(Mobj, 'read_obc_elems')
    oElems = [Mobj.read_obc_elems{:}];
else
    oElems = [];
end
oNodes = [Mobj.read_obc_nodes{:}];
for obc_i = 1:length(oNodes)
    % Find the 20 closest ones.
    [~, h_idx] = sort(sqrt((lon(:) - Mobj.lon(oNodes(obc_i))).^2 + ...
        (lat(:) - Mobj.lat(oNodes(obc_i))).^2));
    obc_hycom_idx = [obc_hycom_idx; h_idx(1:20)];
end
for obc_i = 1:length(oElems)
    % Find the 20 closest ones.
    [~, h_idx] = sort(sqrt((lon(:) - Mobj.lonc(oElems(obc_i))).^2 + ...
        (lat(:) - Mobj.latc(oElems(obc_i))).^2));
    obc_hycom_idx = [obc_hycom_idx; h_idx(1:20)];
end
obc_hycom_idx = unique(obc_hycom_idx, 'stable');
[xidx, yidx] = ind2sub(size(lon), obc_hycom_idx);

% Calculate resolution of parent coarse model
hdx = max([max(diff(lon(:,1))),max(diff(lat(1,:)))]);

% Number of sigma layers.
fz = size(Mobj.siglayz, 2);

% Make a 3D array of the HYCOM depths and mask where we don't have data.
% This can then be used in the interpolation instead of trying to deal with
% this on the fly.
hdepth = permute(repmat(hycom.Depth.data, [1, nx, ny]), [2, 3, 1]);
mask = hycom.(fields{ff}).data(:, :, :, 1) > 1.26e29;
% Used to use the landmask to check whether to extrapolate data onto land.
% Realised this is no longer necessary, so comment this out for future
% removal.
% hdepth(mask) = nan;
% landmask = hycom.(fields{ff}).data(:, :, 1, 1) > 1.26e29;

% Ignore depth layers in the HYCOM data which are below the maximum depth
% of the open boundary nodes so we don't waste time trying to interpolate
% those data. Add a buffer of one layer in case I'm an idiot.
if ~isfield(Mobj, 'hc')
    Mobj.hc = nodes2elems(Mobj.h, Mobj);
end
max_obc_depth = max([max(Mobj.h([Mobj.read_obc_nodes{:}])), max(Mobj.hc([Mobj.read_obc_elems{:}]))]);
[~, zidx] = min(abs(hycom.Depth.data - max_obc_depth));
zidx = min([zidx + 1, length(hycom.Depth.data)]);

if ftbverbose
    tic
end
% Only do warnings about removing values outside some ranges once per
% variable.
warned = true(4, 1);
for v = 1:length(fields)

    if ~(isfield(hycom.(fields{v}), 'data') && ndims(hycom.(fields{v}).data) > 3)
        continue
    end

    % Supply FVCOM grid positions depending on whether we're working on the
    % velocity data (on elements) or anything else (on nodes).
    if any(strcmpi(fields{v}, {'u', 'v'}))
        fvlon = Mobj.lonc(oElems);
        fvlat = Mobj.latc(oElems);
        sigma = Mobj.siglayzc;
        grid_ids = oElems;
        % Number of boundary elements to process
        nf = length(oElems);

        if ftbverbose
            fprintf('Variable %s on elements (%d positions)\n', fields{v}, length(fvlon))
        end
    else
        % Make sure the nodes are listed in the same way as in
        % casename_obc.dat.
        fvlon = Mobj.lon(oNodes);
        fvlat = Mobj.lat(oNodes);
        sigma = Mobj.siglayz;
        grid_ids = oNodes;
        % Number of boundary nodes
        nf = length(oNodes);

        if ftbverbose
            fprintf('Variable %s on nodes (%d positions)\n', fields{v}, length(fvlon))
        end
    end

    fvtemp = nan(nf, fz, nt); % FVCOM interpolated data

    for t = 1:nt

        if ftbverbose
            fprintf('%s : %i of %i %s timesteps... \n', subname, t, nt, fields{v})
        end
        % Get the current 3D array of HYCOM results. Only load the data we
        % need in both the vertical and horizontal (i.e. clip to the
        % locations which cover the open boundary positions only).
        pctemp3 = nan(length(xidx), zidx, 1);
        for oidx = 1:length(xidx)
            pctemp3(oidx, :) = hycom.(fields{v}).data(xidx(oidx), yidx(oidx), 1:zidx, t);
        end
        % pcolor(lon(xidx(1):xidx(2), yidx(1):yidx(2)), lat(xidx(1):xidx(2), yidx(1):yidx(2)), pctemp3(:, :, 1)); shading flat; axis('equal', 'tight'); colorbar
        % hold on
        % plot(fvlon, fvlat, 'ro')
        % Remove NaNs since we've assumed later on that everything above
        % some value is NaN instead (currently 1.26e29). Make the value
        % that we replace the NaNs with larger since we check for larger
        % values when deleting invalid values.
        pctemp3(isnan(pctemp3)) = 1.27e29;

        % Preallocate the intermediate results array.
        itempz = nan(nf, nz);

        for j = 1:nz
            if j == zidx + 1
                if ftbverbose
                    fprintf('   skipping remaining coarse data layers (%d) as their depths are below the FVCOM open boundary depths\n', nz - zidx)
                end
                continue
            elseif j > zidx
                continue
            else
                if ftbverbose
                    fprintf('   coarse data layer %d of %d\n', j, nz)
                end
            end

            % Now extract the relevant layer from the 3D subsets.
            tpctemp2 = pctemp3(:, j);

            % Create and apply a mask to remove values outside the domain.
            % This inevitably flattens the arrays, but it shouldn't be a
            % problem since we'll be searching for the closest values in
            % such a manner as is appropriate for an unstructured grid
            % (i.e. we're assuming the HYCOM data is irregularly sampled
            % and interpolating with a triangulation).
            mask = tpctemp2 > 1.26e29;

            % We need to do some more checks for the data which has been
            % saved via Python. This is a sort of bounds check to eliminate
            % unrealistic data. Warn if we actually delete anything.
            switch fields{v}
                case 'salinity'
                    mask_alt = tpctemp2 < 0;
                    if min(mask_alt(:)) == 1 && warned(1)
                        warned(1) = false;
                        warning('Removing negative salinities from the HYCOM data.')
                    end
                case 'temperature'
                    mask_alt = tpctemp2 < -20;
                    if min(mask_alt(:)) == 1 && warned(2)
                        warned(2) = false;
                        warning('Removing temperature values below -20 celsius from the HYCOM data.')
                    end
                case 'ssh'
                    mask_alt = tpctemp2 < -20;
                    if min(mask_alt(:)) == 1 && warned(3)
                        warned(3) = false;
                        warning('Removing sea surface height values below -20m from the HYCOM data.')
                    end
                case {'u', 'v'}
                    mask_alt = tpctemp2 > 100 | tpctemp2 < -100 ;
                    if min(mask_alt(:)) == 1 && warned(4)
                        warned(3) = false;
                        warning('Removing non-tidal velocities above/below +/-100m/s from the HYCOM data.')
                    end
                otherwise
                    % Some other variable we won't mask.
                    mask_alt = true(size(tpctemp2));
            end
            mask = logical(~(~mask .* ~mask_alt));
            clear mask_alt
            tpctemp2(mask) = [];

            % Also apply the masks to the position arrays so we can't even
            % find positions outside the domain, effectively meaning if a
            % value is outside the domain, the nearest value to the current
            % boundary position will be used.
            % Extract the longitude and latitude data for the open boundary
            % region.
            tlon = lon(obc_hycom_idx);
            tlat = lat(obc_hycom_idx);
            tlon(mask) = [];
            tlat(mask) = [];

            % If the HYCOM depths are now below the maximum depth in the
            % model domain, we'll have no valid data, so skip the
            % interpolation and just leave the NaNs in the itempz array.
            if isempty(tlon) && isempty(tlat)
                continue
            end

            % Preallocate the intermediate results array.
            itempobc = nan(nf, 1);

            % Speed up the tightest loop with a parallelized loop.
            parfor i = 1:nf
                fx = fvlon(i);
                fy = fvlat(i);

                [dist, ii] = sort(sqrt((tlon - fx).^2 + (tlat - fy).^2));
                % Get the n nearest positions from HYCOM data (more?
                % fewer?).
                np = 16;
                if length(ii) < np
                    % Reset np to the number of points we actually have.
                    np = length(ii);
                end
                ixy = ii(1:np);

                % If the minimum distance away is greater than three times
                % the HYCOM grid resolution, skip this point at this
                % vertical level. If this is not done, we can be grabbing
                % data from a significant distance away
                if min(dist) > 3 * hdx
                    continue
                end

                % Get the variables into static variables for the
                % parallelisation.
                plon = double(tlon(ixy));
                plat = double(tlat(ixy));
                ptemp = tpctemp2(ixy);

                % Use a triangulation to do the horizontal interpolation if
                % we have enough points, otherwise take the mean of the two
                % values.
                if length(plon) < 3
                    itempobc(i) = mean(ptemp);
                else
                    try
                        tritemp = scatteredInterpolant(plon, plat, ptemp, 'natural', 'none');
                        itempobc(i) = tritemp(fx, fy);
                    catch err
                        if strcmp(err.identifier, 'MATLAB:subsassignnumelmismatch')
                            warning(['Scatter points failed the', ...
                                ' triangular interpolation. Falling', ...
                                ' back to inverse distance weighted', ...
                                ' interpolation.'])
                            % Use the inverse distance weighted mean of the
                            % values for the interpolated value (the values
                            % are in a straight line and can't be
                            % interpolated with a triangular
                            % interpolation).
                            w = 1 ./ sqrt((plon - fx).^2 + (plat - fy).^2);
                            w = w / max(w);
                            itempobc(i) = sum(ptemp .* w) ./ sum(w);
                        else
                            error(err.message)
                        end
                    end
                end

                % Check we're interpolating at the surface correctly.
                % figure(100)
                % if i == 1
                %     clf
                %     pcolor(hycom.lon - (hdx/2), hycom.lat - (hdx/2), hycom.(fields{v}).data(:, :, j, t)); colorbar; hold on
                %     shading flat
                %     hold on
                % end
                % plot(tlon(ixy), tlat(ixy), 'ko')
                % plot(fx, fy, 'rx')
                % plot(tlon(ixy(1)), tlat(ixy(1)), 'gs')
                % scatter(tlon(ixy(1)), tlat(ixy(1)), 20, ptemp(1), 'filled', 'markeredgecolor', 'r')
                % scatter(fx, fy, 20, ptemp(1), 'filled', 'markeredgecolor', 'g')
                % disp(ptemp(1))
                % caxis([-0.1, 0.1])
                % axis('equal', 'tight')
                % axis([fx - 0.5, fx + 0.5, fy - 0.5, fy + 0.5])
                % legend('Data', 'Candidate valid', 'FVCOM node', 'Nearest valid', 'Interpolated', 'Location', 'NorthOutside', 'Orientation', 'Horizontal')
                % legend('BoxOff')
                % pause(0.1)

                if isnan(itempobc(i))
                    % Use the surface layer as the canonical land mask and
                    % check that the issue here is not just that the
                    % current open boundary position is shallower than this
                    % layer's depth. In the case where we're at the
                    % surface, we always want a value, so use the closest
                    % value, otherwise we can skip this data and leave it
                    % as NaN. The vertical interpolation will strip out the
                    % NaN depths so we shouldn't have any problems from
                    % that perspective.

                    % I used to check we were on land, but really,
                    % itempobc(i) will only equal NaN if we're on land for
                    % this layer. This is only a problem when we're at the
                    % surface as we always need at least one value to do
                    % the vertical interpolation. So, check if we're at the
                    % surface and if so, grab the nearest point. Otherwise,
                    % leave the NaN in place.
                    %[~, jj] = min(sqrt((lon(:) - fx).^2 + (lat(:) - fy).^2));
                    %[ir, ic] = ind2sub(size(lon), jj);
                    %if landmask(ir, ic) == 1 && j == 1
                    if j == 1
                        %fprintf('Sea surface or on land (j = %i, lon: %.5f, lat: %.5f)\n', j, lon(ir, ic), lat(ir, ic))
                        itempobc(i) = tpctemp2(ii(1));

                        %clf
                        %pcolor(hycom.lon, hycom.lat, landmask * 1); colorbar; hold on
                        %plot(lon(ir, ic), lat(ir, ic), 'ko')
                        %plot(fx, fy, 'rx')
                        %plot(tlon(ii(1)), tlat(ii(1)), 'gs')
                        %axis('square')
                        %axis([fx - 1.5, fx + 1.5, fy - 1.5, fy + 1.5])
                        %legend('Land mask', 'Mask test', 'FVCOM node', 'Nearest valid', 'Location', 'NorthOutside', 'Orientation', 'Horizontal')
                        %legend('BoxOff')
                    end
                end
            end

            % Put the results in the intermediate array.
            itempz(:, j) = itempobc;

        end

        % If you want to check the interpolation has worked properly:
        % xx = repmat(fvlon, [1, nz]);
        % yy = repmat(fvlat, [1, nz]);
        % zz = repmat(-hycom.Depth.data, [1, nf])';
        % dd = nanmax(hdepth, [], 3);
        % cc = itempz;
        % mm = ~isnan(cc);
        % ffx = repmat(fvlon, [1, fz]);
        % ffy = repmat(fvlat, [1, fz]);
        % ff = Mobj.siglayz(oNodes, :);
        % figure(10)
        % clf
        % scatter3(xx(mm), yy(mm), zz(mm), 40, cc(mm), 'filled')
        % hold on
        % scatter3(hycom.lon(:), hycom.lat(:), -dd(:), 40, 'c.')
        % scatter3(ffx(:), ffy(:), ff(:), 'k.')
        % colorbar
        % zlim([-300, 0])

        % Now we've interpolated in space, we can interpolate the z-values
        % to the sigma depths.

        fvtempz = nan(nf, fz);
        for pp = 1:nf
            % Get the FVCOM depths at this position
            tfz = sigma(grid_ids(pp), :);  % grid_ids is oNodes or oElems

            % The HYCOM data is unusual in that the depths are fixed and
            % data are only saved at the depths shallower than the
            % bathymetry. As such, we get NaN values below the water depth
            % and we need to filter those out here.

            % Find the HYCOM depths which cover the modelled depth range.
            tpz = -hycom.Depth.data;
            % Mask the HYCOM depths with the data array at this position.
            mm = isnan(itempz(pp, :));
            tpz(mm) = [];

            % If HYCOM has a single value, just repeat it across all depth
            % values.
            if length(tpz) == 1
                fvtempz(pp, :) = repmat(itempz(pp, ~mm), [1, length(tfz)]);
            else
                % To ensure we get the full vertical expression of the
                % vertical profiles, we need to normalise the HYCOM and
                % FVCOM depths to the same range. This is because in
                % instances where FVCOM depths are shallower (e.g. in
                % coastal regions), if we don't normalise the depths, we
                % end up truncating the vertical profile. This approach
                % ensures we always use the full vertical profile, but
                % we're potentially squeezing it into a smaller depth.
                A = max(tpz);
                B = min(tpz);
                C = max(tfz);
                D = min(tfz);
                norm_tpz = (((D - C) * (tpz - A)) / (B - A)) + C;

                % Get the temperature and salinity values for this position
                % and interpolate down the water column (from HYCOM to
                % FVCOM).
                if any(~isnan(norm_tpz))
                    fvtempz(pp, :) = interp1(norm_tpz, itempz(pp, ~mm), tfz, 'pchip', 'extrap');

                    %figure(800);
                    %clf
                    %plot(itempz(pp, ~mm), tpz, 'r-o')
                    %hold on
                    %plot(fvtempz(pp, :), tfz, 'k-x')
                    %legend('HYCOM', 'FVCOM')
                else
                    warning('Should never see this... ') % because we test for NaNs when fetching the values.
                    warning('FVCOM boundary position at %f, %f is outside the HYCOM domain. Skipping.', fvlon(pp), fvlat(pp))
                    continue
                end
            end
        end

        % Find and remove NaNs.
        parfor pp = 1:fz
            test = fvtempz(:, pp);
            if any(isnan(test))
                igood = ~isnan(test);
                ftri = scatteredInterpolant(fvlon(igood), fvlat(igood), test(igood), 'nearest', 'nearest');
                fvtempz(:, pp) = ftri(fvlon,fvlat);
            end
        end

        % The horizontally- and vertically-interpolated values in the final
        % FVCOM results array.
        fvtemp(:, :, t) = fvtempz;

        if ftbverbose
            fprintf('done.\n')
        end
    end
    % Dump the data into a temporary structure.
    fvcom.(fields{v}) = fvtemp;
end % loop through variables
if ftbverbose
    toc
end

fvfields = fieldnames(fvcom);
for s = 1:length(fvfields)
    switch fvfields{s}
        case 'temperature'
            Mobj.temperature = fvcom.temperature;
        case 'salinity'
            Mobj.salt = fvcom.salinity;
        case 'u'
            Mobj.u = fvcom.u;
        case 'v'
            Mobj.v = fvcom.v;
        case 'density'
            Mobj.rho1 = fvcom.density;
        case 'ssh'
            Mobj.ssh = fvcom.ssh;
        otherwise
            warning('Unrecognised variable %s', fvfields{s})
    end
end

if isfield(hycom, 'time')
    Mobj.ts_times = hycom.time;
end

cleaner = onCleanup(@() delete(gcp('nocreate')));

if ftbverbose
    fprintf('end   : %s\n', subname)
end


%%
% Plot a vertical profile for a boundary node (for my Irish Sea case, this
% is one of the ones along the Celtic Sea boundary). Also plot the
% distribution of interpolated values over the HYCOM data. Add the location
% of the vertical profile (both FVCOM and HYCOM) to the plot.

% nn = 128;  % open boundary index
% tt = 1;    % time index
% fvz = 1;   % fvcom depth index (currently 1-20)
% hyz = 1;   % hycom depth index (1-33)
%
% % Find the HYCOM seabed indices
% % [~, hyz] = nanmax(hdepth, [], 3);
%
% % Get the corresponding indices for the HYCOM data
% [~, idx] = min(sqrt((lon(:) - fvlon(nn)).^2 + (lat(:) - fvlat(nn)).^2));
% [xidx, yidx] = ind2sub(size(lon), idx);
%
% zidx = isfinite(hdepth(xidx, yidx, :));
% hz = 1:nz;
%
% % close all
%
% figure(100)
% clf
% subplot(2,2,1)
% plot(Mobj.temperature(nn, :, tt), Mobj.siglayz(oNodes(nn), :), 'x-')
% xlabel('Temperature (^{\circ}C)')
% ylabel('Depth (m)')
% title('FVCOM')
%
% subplot(2,2,2)
% plot(squeeze(hycom.temperature.data(xidx, yidx, zidx, tt)), squeeze(-hdepth(xidx, yidx, zidx)), 'rx-')
% xlabel('Temperature (^{\circ}C)')
% ylabel('Depth (m)')
% title('HYCOM')
%
% subplot(2,2,3)
% plot(Mobj.temperature(nn, :, tt), 1:fz, 'x-')
% xlabel('Temperature (^{\circ}C)')
% ylabel('Array index')
% title('FVCOM')
%
% subplot(2,2,4)
% plot(squeeze(hycom.temperature.data(xidx, yidx, zidx, tt)), hz(zidx), 'rx-')
% xlabel('Temperature (^{\circ}C)')
% ylabel('Array index')
% title('HYCOM')
%
% % Figure to check everything's as we'd expect. Plot first time step with
% % the HYCOM surface temperature as a background with the interpolated
% % boundary node surface values on top.
%
% figure(200)
% clf
% % Plot HYCOM surface data (last sigma layer)
% dx = mean(diff(hycom.lon(:)));
% dy = mean(diff(hycom.lat(:)));
% temp = hycom.temperature.data(:, :, :, tt);
% pcolor(hycom.lon - (dx / 2), hycom.lat - (dy / 2), ...
%     squeeze(temp(:, :, hyz)))
% shading flat
% axis('equal', 'tight')
% daspect([1.5, 1, 1])
% hold on
% % Add the interpolated surface data (first sigma layer)
% scatter(Mobj.lon(oNodes), Mobj.lat(oNodes), 40, Mobj.temperature(:, fvz, tt), 'filled', 'MarkerEdgeColor', 'k')
% axis([min(Mobj.lon(oNodes)), max(Mobj.lon(oNodes)), min(Mobj.lat(oNodes)), max(Mobj.lat(oNodes))])
% caxis([3, 13])
% plot(lon(xidx, yidx), lat(xidx, yidx), 'rs') % polcoms is all backwards
% plot(Mobj.lon(oNodes(nn)), Mobj.lat(oNodes(nn)), 'wo')
% colorbar
