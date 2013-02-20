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
%               - Mobj.obc_nodes - list of open boundary node inidices.
%               - Mobj.nObcNodes - number of nodes in each open boundary.
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
%    Mobj = MATLAB structure in which two new fields (called meanflow,
%           and mf_times). meanflow has a size of (sum(Mobj.nObcNodes),
%           time). The time dimension is determined based on the input
%           NetCDF file. The mf_time variable is just the input file times
%           in Modified Julian Day.
%
% EXAMPLE USAGE
%    Mobj = get_POLCOMS_meanflow(Mobj, files)
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2013-02-20 First version.
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

fvlon = Mobj.lon(Mobj.obc_nodes(Mobj.obc_nodes ~= 0));
fvlat = Mobj.lat(Mobj.obc_nodes(Mobj.obc_nodes ~= 0));

% Number of boundary nodes
nf = sum(Mobj.nObcNodes);

fvmf = nan(nf, nt); % FVCOM interpolated mean flow

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

    % Calculate the velocity and direction (?) from the u and v vectors.
    % This approach means we don't have to interpolate u and v separately
    % (which is expensive) but just do velocity on its own.
    pcvel3 = sqrt(pcu3.^2 + pcv3.^2);
       
    % For the velocity, it appears we don't need to have a depth varying
    % value (instead mean flow is scaled according to the values in MFDIST,
    % although I'm not yet sure exactly what those values are...). So,
    % we'll calculate a depth average velocity here for the time being.
    % Transpose the array here so it's (x, y) rather than the original (y,
    % x).
    pcvel3mean = mean(pcvel3, 3)';

    % We need to create a mask to eliminate land values and apply it to the
    % depth averaged values.
    mask = squeeze(pc.depth.data(:, :, end, t))' >= 0;
    
    tpcvel3mean = pcvel3mean;
    tlon = lon;
    tlat = lat;
    % Clear out the land values and flatten the arrays.
    tpcvel3mean(mask) = [];
    tlon(mask) = [];
    tlat(mask) = [];
    
    % Speed up the tightest loop with a parallelized loop.
    parfor i = 1:nf
        % Now we can do each position within the depth averaged array.

        fx = fvlon(i);
        fy = fvlat(i);

        [~, ii] = sort(sqrt((tlon - fx).^2 + (tlat - fy).^2));
        % Get the n nearest nodes from PML POLCOMS-ERSEM data (more?
        % fewer?).
        ixy = ii(1:16);

        % Get the variables into static variables for the
        % parallelisation.
        plon = tlon(ixy);
        plat = tlat(ixy);
        pvel = tpcvel3mean(ixy);

        % Use a triangulation to do the horizontal interpolation.
        trivel = TriScatteredInterp(plon', plat', pvel', 'natural');
        ivelobc(i) = trivel(fx, fy);

        % Check if we have NaNs (mostly if the position is outside the
        % model domain).
        if isnan(ivelobc(i))
            warning('FVCOM boundary node at %f, %f is outside the PML POLCOMS-ERSEM domain. Setting to the closest PML POLCOMS-ERSEM value.', fx, fy)
            ivelobc(i) = tpcvel2(ii(1));
        end
    end

  
%     for j = 1:nz
%         % Now extract the relevant layer from the 3D subsets. Transpose the
%         % data to be (x, y) rather than (y, x).
%         pcvel2 = pcvel3(:, :, j)';
%         pcdepth2 = squeeze(pc.depth.data(:, :, j, t))';
% 
%         % Create new arrays which will be flattened when masking (below).
%         tpcvel2 = pcvel2;
%         tpcdepth2 = pcdepth2;
%         tlon = lon;
%         tlat = lat;
%         
%         % Create and apply a mask to remove values outside the domain. This
%         % inevitably flattens the arrays, but it shouldn't be a problem
%         % since we'll be searching for the closest values in such a manner
%         % as is appropriate for an unstructured grid (i.e. we're assuming
%         % the PML POLCOMS-ERSEM data is irregularly spaced).
%         mask = tpcdepth2 > 20000;
%         tpcvel2(mask) = [];
%         tpcdepth2(mask) = [];
%         % Also apply the masks to the position arrays so we can't even find
%         % positions outside the domain, effectively meaning if a value is
%         % outside the domain, the nearest value to the boundary node will
%         % be used.
%         tlon(mask) = [];
%         tlat(mask) = [];
% 
%         % Preallocate the intermediate results arrays.
%         ivelobc = nan(nf, 1);
%         idepthobc = nan(nf, 1);
%         
%         % Speed up the tightest loop with a parallelized loop.
%         parfor i = 1:nf
%             % Now we can do each position within the 2D layer.
% 
%             fx = fvlon(i);
%             fy = fvlat(i);
% 
%             [~, ii] = sort(sqrt((tlon - fx).^2 + (tlat - fy).^2));
%             % Get the n nearest nodes from PML POLCOMS-ERSEM data (more?
%             % fewer?).
%             ixy = ii(1:16);
% 
%             % Get the variables into static variables for the
%             % parallelisation.
%             plon = tlon(ixy);
%             plat = tlat(ixy);
%             pvel = tpcvel2(ixy);
%             pdepth = tpcdepth2(ixy);
%             
%             % Use a triangulation to do the horizontal interpolation.
%             trivel = TriScatteredInterp(plon', plat', pvel', 'natural');
%             triz = TriScatteredInterp(plon', plat', pdepth', 'natural');
%             ivelobc(i) = trivel(fx, fy);
%             idepthobc(i) = triz(fx, fy);
%             
%             % Check both, though if one is NaN, they both will be.
%             if isnan(ivelobc(i)) || isnan(idepthobc(i))
%                 warning('FVCOM boundary node at %f, %f is outside the PML POLCOMS-ERSEM domain. Setting to the closest PML POLCOMS-ERSEM value.', fx, fy)
%                 ivelobc(i) = tpcvel2(ii(1));
%                 idepthobc(i) = tpcdepth2(ii(1));
%             end
%         end
%         
%         % Put the results in this intermediate array.
%         ivelz(:, j) = ivelobc;
%         idepthz(:, j) = idepthobc;
%     end

    % The horizontally-interpolated values in the final FVCOM results
    % array.
    fvmf(:, t) = ivelobc;

    if ftbverbose
        fprintf('done.\n')
    end
end
if ftbverbose
    toc
end

% Stick the values in the mesh structure.
Mobj.velocity = fvmf;

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

        
