function Mobj = get_POLCOMS_meanflow(Mobj, files)
% Read mean flow from the PML POLCOMS-ERSEM NetCDF AMM model output files
% and interpolate onto the open boundaries in Mobj.
%
% function Mobj = get_POLCOMS_meanflow(Mobj, ts, polcoms_bathy, varlist)
%
% DESCRIPTION:
%    Interpolate u and v flow vectors to calculate daily mean flow at the
%    FVCOM open boundaries at all sigma levels.
%
% INPUT:
%   Mobj    = MATLAB mesh structure which must contain:
%               - Mobj.lon, Mobj.lat - node coordinates (lat/long).
%               - Mobj.obc_elements - list of open boundary element
%               inidices.
%               - Mobj.nObcElements - number of elements in each open boundary.
%   files   = Cell array of PML POLCOMS-ERSEM NetCDF file(s) in which 4D
%             variables of u and v velocity components (called 'ucurD' and
%             'vcurD') exist. Their shape should be (y, x, sigma, time).
%
% NOTES:
%
%   - If you supply multiple files in files, there are a few assumptions:
%
%       - Variables are only appended if there are 4 dimensions; fewer than
%       that, and the values are assumed to be static across all the given
%       files (e.g. longitude, latitude etc.). The fourth dimension is
%       time.
%       - The order of the files given should be chronological.
% 
%   - The NetCDF files used here are those from the PML POLCOMS-ERSEM model
%   output.
%
% OUTPUT:
%    Mobj = MATLAB structure in which six new fields (mf_time, meanflow_u,
%           meanflow_v, meanflow_ubar, meanflow_ubar and velocity).
%           meanflow_ubar, meanflow_vbar and velocity have sizes of
%           (sum(Mobj.nObcNodes), time); meanflow_u and meanflow_v have
%           sizes (sum(Mobj.nObcNodes), siglay, time). The time dimension
%           is determined based on the input NetCDF file. The mf_times
%           variable is just the input file times in Modified Julian Day.
%
% EXAMPLE USAGE
%    Mobj = get_POLCOMS_meanflow(Mobj, files)
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2013-02-20 First version.
%    2013-02-26 Add interpolation of the u and v vectors separately and
%    then calculate the interpolated velocity at the end.
%    2013-02-27 Change the vertical interpolation to be scaled within the
%    POLCOMS-ERSEM depth range for the current node. The net result is that
%    the vertical profiles are squashed or stretched to fit within the
%    FVCOM depths. This means the full profile structure is maintained in
%    the resulting FVCOM boundary input despite the differing depths at the
%    FVCOM boundary node.
%    2013-02-28 Change the interpolation to occur on the open boundary
%    elements rather than on the open boundary nodes. This requires a
%    couple of extra steps in the pre-processing (notably running
%    find_boundary_elements) as well as transferring the sigma depths to
%    the element centres.
%
%==========================================================================

subname = 'get_POLCOMS_meanflow';

global ftbverbose;
if ftbverbose
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end

varlist = {'lon', 'lat', 'ucurD', 'vcurD', 'time', 'depth', 'pdepthD'};

pc = get_POLCOMS_netCDF(files, varlist);

[~, ~, nz, nt] = size(pc.ucurD.data);

[lon, lat] = meshgrid(pc.lon.data, pc.lat.data);

% We need positions of the open boundary element centroids (i.e. the MCE
% centroid) at which to calculate the mean flow velocity. As such, we need
% to first calculate a number of parameters at the element centres from the
% nodal values we typically collect (e.g. positions, sigma layers etc.). Do
% all that here.
lonc = nodes2elems(Mobj.lon, Mobj);
latc = nodes2elems(Mobj.lat, Mobj);
% For the sigma levels, we need to wrap the conversion in a loop.
siglayzc = nan([Mobj.nElems, size(Mobj.siglayz, 2)]);
for zz = 1:size(Mobj.siglayz, 2)
    siglayzc(:, zz) = nodes2elems(Mobj.siglayz(:, zz), Mobj);
end

oElements = [Mobj.read_obc_elements{:}];

% Find the FVCOM positions for all the open boundary elements.
fvlon = lonc(oElements);
fvlat = latc(oElements);

% Number of open boundary elements
ne = sum(Mobj.nObcElements);
% Number of sigma layers.
fz = size(siglayzc, 2);

fvmfu = nan(ne, fz, nt); % FVCOM interpolated flow vector components
fvmfv = nan(ne, fz, nt); % FVCOM interpolated flow vector components

if ftbverbose
    tic
end

for t = 1:nt
    if ftbverbose
        fprintf('%s : %i of %i timesteps... ', subname, t, nt)
    end
    % Get the current 3D array of PML POLCOMS-ERSEM results.
    pcu3 = pc.ucurD.data(:, :, :, t);
    pcv3 = pc.vcurD.data(:, :, :, t);
    pcz3 = pc.depth.data(:, :, :, t);

    % Preallocate the intermediate results arrays.
    iuz = nan(ne, nz);
    ivz = nan(ne, nz);
    izz = nan(ne, nz);

    % We need to create a mask to eliminate land values and apply it to the
    % depth averaged values.
    mask = squeeze(pc.depth.data(:, :, end, t))' >= 0;

    % Mask the longitude and latitude here (since it's depth independent,
    % there's no point doing it inside the next loop).
    tlon = lon;
    tlat = lat;
    tlon(mask) = [];
    tlat(mask) = [];

    for z = 1:nz
        % Turns out we do need vertical velocity for the mean flow, so
        % iterate through each vertical layer before interpolating
        % horizontally. The vertical interpolation happens after the
        % horizontal has been completed. Transpose the arrays to be (x, y)
        % rather than (y, x).
        pcu2 = pcu3(:, :, z)';
        pcv2 = pcv3(:, :, z)';
        pcz2 = pcz3(:, :, z)';

        % Flatten the arrays through the masking.
        pcu2(mask) = [];
        pcv2(mask) = [];
        pcz2(mask) = [];

        % Preallocate the intermediate results arrays.
        iuobc = nan(ne, 1);
        ivobc = nan(ne, 1);
        izobc = nan(ne, 1);

        % Speed up the tightest loops with a parallelized loop.
        parfor i = 1:ne
            % Now we can do each position within the current depth layer.

            fx = fvlon(i);
            fy = fvlat(i);

            [~, ii] = sort(sqrt((tlon - fx).^2 + (tlat - fy).^2));
            % Get the n nearest elements from PML POLCOMS-ERSEM data (more?
            % fewer?).
            ixy = ii(1:16);

            % Get the variables into static variables for the
            % parallelisation.
            plon = tlon(ixy);
            plat = tlat(ixy);
            pu = pcu2(ixy);
            pv = pcv2(ixy);
            pz = pcz2(ixy);

            % Use a triangulation to do the horizontal interpolation.
            triu = TriScatteredInterp(plon', plat', pu', 'natural');
            triv = TriScatteredInterp(plon', plat', pv', 'natural');
            triz = TriScatteredInterp(plon', plat', pz', 'natural');
            iuobc(i) = triu(fx, fy);
            ivobc(i) = triv(fx, fy);
            izobc(i) = triz(fx, fy);

            % Check if we have NaNs (mostly if the position is outside the
            % model domain).
            if isnan(iuobc(i)) || isnan(ivobc(i))
                warning('FVCOM boundary element at %f, %f is outside the PML POLCOMS-ERSEM domain. Setting to the closest PML POLCOMS-ERSEM value.', fx, fy)
                iuobc(i) = pcu2(ii(1));
                ivobc(i) = pcv2(ii(1));
                izobc(i) = pcz2(ii(1));
            end
        end

        % Put the results in this intermediate array.
        iuz(:, z) = iuobc;
        ivz(:, z) = ivobc;
        izz(:, z) = izobc;
    end

    % Now we've interpolated in space, we can interpolate the z-values
    % to the sigma depths.

    % Preallocate the output arrays
    fvuz = nan(ne, fz);
    fvvz = nan(ne, fz);

    for pp = 1:ne
        % Get the FVCOM depths at this element
        tfz = siglayzc(oElements(pp), :);
        % Now get the interpolated PML POLCOMS-ERSEM depth at this element
        tpz = izz(pp, :);

        % To ensure we get the full vertical expression of the vertical
        % profiles, we need to normalise the POLCOMS-ERSEM depths to the
        % FVCOM range for the current element. This is because in instances
        % where FVCOM depths are shallower (e.g. in coastal regions), if we
        % don't normalise the depths, we end up truncating the vertical
        % profile. This approach ensures we always use the full vertical
        % profile, but we're potentially squeezing it into a smaller depth.
        A = max(tpz);
        B = min(tpz);
        C = max(tfz);
        D = min(tfz);
        norm_tpz = (((D - C) * (tpz - A)) / (B - A)) + C;

        % Get the u and v velocity values for this elment and interpolate
        % down the water column (from PML POLCOMS-ERSEM to FVCOM). I
        % had originally planned to use csaps for the vertical
        % interplation/subsampling at each location. However, the demo
        % of csaps in the MATLAB documentation makes the interpolation
        % look horrible (shaving off extremes). interp1 provides a
        % range of interpolation schemes, of which pchip seems to do a
        % decent job of the interpolation (at least qualitatively).
        if ~isnan(tpz)
            fvuz(pp, :) = interp1(norm_tpz, iuz(pp, :), tfz, 'linear', 'extrap');
            fvvz(pp, :) = interp1(norm_tpz, ivz(pp, :), tfz, 'linear', 'extrap');
        else
            warning('Should never see this... ') % because we test for NaNs when fetching the values.
            warning('FVCOM boundary element at %f, %f is outside the PML POLCOMS-ERSEM domain. Skipping.', fvlon(pp), fvlat(pp))
            continue
        end
    end

    % The horizontally- and vertically-interpolated values in the final
    % FVCOM results array.
    fvmfu(:, :, t) = fvuz;
    fvmfv(:, :, t) = fvvz;

    if ftbverbose
        fprintf('done.\n')
    end
end
if ftbverbose
    toc
end

% Stick the values in the mesh structure.
Mobj.meanflow_u = fvmfu;
Mobj.meanflow_v = fvmfv;

% Now we have the 3D arrays, create depth averaged velocities too
Mobj.meanflow_ubar = squeeze(mean(Mobj.meanflow_u, 2));
Mobj.meanflow_vbar = squeeze(mean(Mobj.meanflow_v, 2));

% Depth averaged velocity
Mobj.velocity = squeeze(mean(sqrt(Mobj.meanflow_u.^2 + Mobj.meanflow_v.^2), 2));

% Convert the current times to Modified Julian Day (this is a bit ugly).
pc.time.all = strtrim(regexp(pc.time.units, 'since', 'split'));
pc.time.datetime = strtrim(regexp(pc.time.all{end}, ' ', 'split'));
pc.time.ymd = str2double(strtrim(regexp(pc.time.datetime{1}, '-', 'split')));
pc.time.hms = str2double(strtrim(regexp(pc.time.datetime{2}, ':', 'split')));

Mobj.mf_times = greg2mjulian(...
    pc.time.ymd(1), ...
    pc.time.ymd(2), ...
    pc.time.ymd(3), ...
    pc.time.hms(1), ...
    pc.time.hms(2), ...
    pc.time.hms(3)) + (pc.time.data / 3600 / 24);

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end

% Check the interpolation along the boundary at the last time step. N.B.
% Since the depths are scaled from the POLCOMS range to the FVCOM range,
% the profiles won't match perfectly in these plots. This is because the
% interpolated FVCOM profiles use the full water column from the POLCOMS
% data rather than truncating it at the FVCOM depth.

% % Map the open boundaries with the depth averaged velocity as a background
% figure(1)
% clf
% imagesc(pc.lon.data, pc.lat.data, mean(sqrt(pcu3.^2 + pcv3.^2), 3)')
% shading flat
% axis('equal', 'tight')
% set(gca,'YDir','normal')
% colorbar
% caxis([0, 0.05])
% hold on
% plot(fvlon, fvlat, 'wo')
% axis([min(fvlon) - 0.1, max(fvlon) + 0.1, min(fvlat) - 0.1, max(fvlat) + 0.1])
%
% for i = 1:ne
%
%     % Add the current element's position and value as a scatter point to
%     % the map plot.
%     scatter(fvlon(i), fvlat(i), 50, Mobj.velocity(i, end), 'filled')
%
%     % Do vertical profiles of u, v and velocity.
%     figure(2)
%     clf
%
%     subplot(3,1,1)
%     % FVCOM vertical profile. Has to be (i, :, end) because the
%     % corresponding POLCOMS data isn't stored as a function of time (i.e.
%     % iuz, ivz and izz are all for the last time step only).
%     fvuz = Mobj.meanflow_u(i, :, end);
%     fvz = siglayzc(oElements(i), :);
%     plot(fvuz, fvz)
%     hold on
%     % The interpolated POLCOMS vertical profile (last time step only)
%     plot(iuz(i, :), izz(i, :), 'g')
%     % The depth-averaged velocity (again, last time step only)
%     fvubar = Mobj.meanflow_ubar(i, end);
%     plot([fvubar, fvubar], [min(fvz), max(fvz)], 'k')
%     xlim([-0.1, 0.2])
%     ylim([-100, 0])
%     xlabel('Mean flow u-velocity (m^{3}s^{-1})')
%     ylabel('Depth (m)')
%     title('u-component')
%     legend('FVCOM', 'POLCOMS', 'Mean FVCOM', 'Location', 'SouthEast')
%     legend('boxoff')
%
%     subplot(3,1,2)
%     % FVCOM vertical profile. Has to be (i, :, end) because the
%     % corresponding POLCOMS data isn't stored as a function of time (i.e.
%     % iuz, ivz and izz are all for the last time step only).
%     fvvz = Mobj.meanflow_v(i, :, end);
%     fvz = siglayzc(oElements(i), :);
%     plot(fvvz, fvz)
%     hold on
%     % The interpolated POLCOMS vertical profile (last time step only)
%     plot(ivz(i, :), izz(i, :), 'g')
%     % The depth-averaged velocity (again, last time step only)
%     fvvbar = Mobj.meanflow_vbar(i, end);
%     plot([fvvbar, fvvbar], [min(fvz), max(fvz)], 'k')
%     xlim([-0.1, 0.2])
%     ylim([-100, 0])
%     xlabel('Mean flow v-velocity (m^{3}s^{-1})')
%     ylabel('Depth (m)')
%     title('v-component')
%     legend('FVCOM', 'POLCOMS', 'Mean FVCOM', 'Location', 'SouthEast')
%     legend('boxoff')
%
%     subplot(3,1,3)
%     % FVCOM vertical profile. Has to be (i, :, end) because the
%     % corresponding POLCOMS data isn't stored as a function of time (i.e.
%     % iuz, ivz and izz are all for the last time step only).
%     fvvelz = sqrt(Mobj.meanflow_u(i, :, end).^2 + Mobj.meanflow_v(i, :, end).^2);
%     fvz = siglayzc(oElements(i), :);
%     plot(fvvelz, fvz)
%     hold on
%     % The interpolated POLCOMS vertical profile (last time step only)
%     plot(sqrt(iuz(i, :).^2 + ivz(i, :).^2), izz(i, :), 'g')
%     % The depth-averaged velocity (again, last time step only)
%     fvvelbar = Mobj.velocity(i, end);
%     plot([fvvelbar, fvvelbar], [min(fvz), max(fvz)], 'k')
%     xlim([-0.1, 0.2])
%     ylim([-100, 0])
%     xlabel('Mean flow velocity (m^{3}s^{-1})')
%     ylabel('Depth (m)')
%     title('velocity')
%     legend('FVCOM', 'POLCOMS', 'Mean FVCOM', 'Location', 'SouthEast')
%     legend('boxoff')
%     pause(0.1)
% end

