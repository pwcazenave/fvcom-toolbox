function dataz = grid_vert_interp(Mobj, lon, lat, data, depth, mask, varargin)
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
%                   - Mobj.have_lonlat - boolean for spherical coordinate
%                   fields (Mobj.lon, Mobj.lat) presence (true = yes).
%                   - Mobj.siglayz - sigma layer depths for all model
%                   nodes.
%   lon, lat    = Rectangular arrays of longitude and latitude (see
%   meshgrid).
%   data, depth = x by y by z (where z is vertical layers) grids of the
%   data and water depths to be interpolated onto the vertical grid defined
%   by Mobj.siglayz.
%   mask        = logical array of positions outside the regularly gridded
%   domain (e.g. if the regular data contains NaNs or other undefined
%   values, create a logical array of those positions so they can be
%   omitted quickly).
%   'extrapolate' [optional] = keyword-argument pair in which the argument
%   specifies the coordinates to use for the extrapolation (e.g.
%   'extrapolate', [Mobj.lonc, Mobj.latc] to extrapolate onto the element
%   centres). Defaults to the element nodes (i.e. [Mobj.lon, Mobj.lat]).
%   
% OUTPUT:
%   x by y by z array of vertically interpolated values at each regular
%   grid location.
%
% EXAMPLE USAGE
%   Basic usage:
%       grid_vert_interp(Mobj, lon, lat, data, depth, mask)
%
%   Extrapolate using the element centres.
%       grid_vert_interp(Mobj, lon, lat, data, depth, mask, ...
%           'extrapolate', [Mobj.lonc, Mobj.latc])
% 
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2013-02-08 First version.
%   2013-05-16 Add support for parallel for-loops (not mandatory, but
%   enabled if the Parallel Computing Toolbox is available).
%   2014-04-28 Add new argument to allow specifying different coordinates
%   for the extrapolation. This allows us to interpolate data which belongs
%   either to the element centres or element nodes (defaults to element
%   nodes). This is only really important when the model grid falls outside
%   the coverage of the supplied data and we're extrapolating data. Update
%   the parallel pool code to use the new parpool function instead of
%   matlabpool in anticipation of the latter's eventual removal from
%   MATLAB. Also update the help.
%   2014-06-12 Fix bug in interpolating in the vertical when HYCOM data has
%   only a single depth bin.
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

% Extract the extrapolation coordinates. Default to nodes but use
% whatever's given if we have the 'extrapolate' argument.
ulon = Mobj.lon;
ulat = Mobj.lat;
for aa = 1:2:length(varargin)
    switch varargin{aa}
        case 'extrapolate'
            ulon = varargin{aa + 1}(:, 1);
            ulat = varargin{aa + 1}(:, 2);
    end
end

wasOpened = false;
if license('test', 'Distrib_Computing_Toolbox')
    % We have the Parallel Computing Toolbox, so launch a bunch of workers.
    try
        % New version for MATLAB 2014a (I think) onwards.
        if isempty(gcp('nocreate'))
            pool = parpool('local');
            wasOpened = true;
        end
    catch
        % Version for pre-2014a MATLAB.
        if matlabpool('size') == 0
            % Force pool to be local in case we have remote pools available.
            matlabpool open local
            wasOpened = true;
        end
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
        [md, mi] = min(sqrt((ulon - lon(xi, yi)).^2 + (ulat - lat(xi, yi)).^2));

        % Skip data point if the closest FVCOM node is more than 10 minutes
        % away.
        if md > 10 / 60
%             if ftbverbose
%                 fprintf('%s : skipping %f, %f (more than 10 minutes from the nearest unstructured grid node),\n', subname, lon(yi, xi), lat(yi, xi))
%             end
            continue
        else
            % Use the FVCOM node's sigma depths to interpolate this regular
            % grid position's temperature and salinity data.
            
            % Get the FVCOM depths closest to this regular grid position.
            try
                tfz = Mobj.siglayz(mi, :);
            catch
                tfz = Mobj.siglayzc(mi, :);
            end
            % Now get the regular grid depths at this node for all the
            % vertical layers and linearly interpolate through the water
            % column onto the FVCOM vertical profile.
            tpz = xdepth(yi, :);

            % Remove any NaN values in the vertical depths (as is the case
            % with the HYCOM data, and interpolate only the data we have
            % that are finite).
            nidx = isnan(tpz);
            % If we have only a single value, repeat it down the entire
            % water column, otherwise, interpolate linearly (with
            % extrapolation).
            if length(tpz(~nidx)) == 1
                ydata(yi, :) = tpz(~nidx);
            else
                ydata(yi, :) = interp1(tpz(~nidx), xdata(yi, ~nidx), tfz, 'linear', 'extrap');
            end
        end
    end
    dataz(xi, :, :) = ydata;
end

% Close the MATLAB pool if we opened it.
if wasOpened
    try
        pool.delete
    catch
        matlabpool close
    end
end

if ftbverbose
    toc
end

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end
