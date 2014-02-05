function dataz = grid_vert_interp(Mobj, lon, lat, data, depth, mask)
% Child function to interpolate a given 3D array (x by y by sigma) values
% to the unstructured grid vertical layers.
%
% function grid_vert_interp(Mobj, lon, lat, data, depth, mask)
%
% DESCRIPTION:
%    Interpolate the regularly gridded data described by lon, lat, data
%    (whose size should be x by y for the lon and lat arrays, and x by y by
%    z for the data array) to the unstructured grid vertical layers defined
%    in Mobj.siglayz. The depths in depth are used to make sure the profile
%    is interpolated the right way up. The closest unstructured grid node
%    is used for the vertical interpolation of each regularly gridded node.
%
% INPUT:
%   Mobj        = MATLAB mesh structure which must contain:
%                   - Mobj.siglayz - sigma layer depths for all model
%                   nodes.
%                   - Mobj.lon, Mobj.lat - node coordinates (long/lat).
%   lon, lat    = Rectangular arrays of longitude and latitude (see
%   meshgrid).
%   data, depth = x by y by z (where z is vertical layers) grids of the
%   data and water depths to be interpolated onto the vertical grid defined
%   by Mobj.siglayz.
%   mask        = logical array of positions outside the regularly gridded
%   domain (e.g. if the regular data contains NaNs or other undefined
%   values, create a logical array of those positions so they can be
%   omitted quickly).
%   
% OUTPUT:
%   x by y by z array of vertically interpolated values at each regular
%   grid location.
%
% EXAMPLE USAGE
%   grid_vert_interp(Mobj, lon, lat, data, depth, mask)
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

subname = 'grid_vert_interp';

global ftbverbose
if ftbverbose
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end

if ~isfield(Mobj, 'siglayz')
    error('Error: missing required sigma layer depth values for the unstructured grid.')
end

if ~Mobj.have_lonlat
    error('Need spherical coordinates')
end

if sum(size(lon) ~= size(data(:, :, 1))) ~= 0 || sum(size(lat) ~= size(data(:, :, 1))) ~= 0
    error('Size of the longitude or latitude arrays do not match the supplied data array')
end

wasOpened = false;
if license('test', 'Distrib_Computing_Toolbox')
    % We have the Parallel Computing Toolbox, so launch a bunch of workers.
    if matlabpool('size') == 0
        % Force pool to be local in case we have remote pools available.
        matlabpool open local
        wasOpened = true;
    end
end

[~, fz] = size(Mobj.siglayz);
[nx, ny, ~, ~] = size(data);

% Preallocate the output arrays
dataz = nan(nx, ny, fz);

if ftbverbose
    tic
end

parfor xi = 1:nx
    % Get all the y and z dimension data for the current x position
    % (temperature, salinity and depth).
    xdata = squeeze(data(xi, :, :));
    xdepth = squeeze(depth(xi, :, :));
    xmask = mask(xi, :);
    
    % Preallocate the arrays for the inner loop.
    ydata = nan(ny, fz);
    for yi = 1:ny
        if xmask(yi)
            continue
        end

        % Find the nearest sigma layer z values from the unstructured grid.
        [md, mi] = min(sqrt((Mobj.lon - lon(xi, yi)).^2 + (Mobj.lat - lat(xi, yi)).^2));

        % Skip data point if the closest FVCOM node is more than 10 minutes
        % away.
        if md > 10 / 60
%             if ftbverbose
%                 fprintf('%s : skipping %f, %f (more than 10 minutes from the nearest unstructured grid node),\n', subname, lon(yi, xi), lat(yi, xi))
%             end
            continue
        else
            % Use the FVCOM node's sigma depths to interpolate this POLCOMS
            % position's temperature and salinity data.
            
            % Get the FVCOM depths closest to this POLCOMS grid position.
            tfz = Mobj.siglayz(mi, :);
            % Now get the POLCOMS depths at this node for all the vertical
            % layers and linearly interpolate through the water column onto
            % the FVCOM vertical profile.
            tpz = xdepth(yi, :);

            % Remove any NaN values in the vertical depths (as is the case
            % with the HYCOM data, and interpolate only the data we have
            % that are finite).
            nidx = isnan(tpz);
            ydata(yi, :) = interp1(tpz(~nidx), xdata(yi, ~nidx), tfz, 'linear', 'extrap');
        end
    end
    dataz(xi, :, :) = ydata;
end

% Close the MATLAB pool if we opened it.
if wasOpened
    matlabpool close
end

if ftbverbose
    toc
end

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end
