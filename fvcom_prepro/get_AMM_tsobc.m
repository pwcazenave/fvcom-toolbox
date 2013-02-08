function Mobj = get_AMM_tsobc(Mobj, ts)
% Read temperature and salinity from AMM NetCDF model output files and
% interpolate onto the open boundaries in Mobj.
%
% function Mobj = get_AMM_tsobc(Mobj, ts, polcoms_bathy, varlist)
%
% DESCRIPTION:
%    Interpolate temperature and salinity values onto the FVCOM open
%    boundaries at all sigma levels.
%
% INPUT:
%   Mobj    = MATLAB mesh structure which must contain:
%               - Mobj.siglayz - sigma layer depths for all model nodes.
%               - Mobj.lon, Mobj.lat - node coordinates (lat/long).
%               - Mobj.obc_nodes - list of open boundary node inidices.
%               - Mobj.nObcNodes - number of nodes in each open boundary.
%   ts      = Cell array of AMM NetCDF file(s) in which 4D
%             variables of temperature and salinity (called 'ETWD' and
%             'x1XD') exist. Their shape should be (y, x, sigma, time).
%
% NOTES:
%
%   - If you supply multiple files in ts, there are a few assumptions:
%
%       - Variables are only appended if there are 4 dimensions; fewer than
%       that, and the values are assumed to be static across all the given
%       files (e.g. longitude, latitude etc.). The fourth dimension is
%       time.
%       - The order of the files given should be chronological.
% 
%   - The NetCDF files used here are those from the MyOcean portal
%   available for free at:
%
%       http://www.myocean.eu/web/24-catalogue.php
%
% OUTPUT:
%    Mobj = MATLAB structure in which three new fields (called temperature,
%           salinity and ts_time). temperature and salinity have sizes
%           (sum(Mobj.nObcNodes), sigma, time). The time dimension is
%           determined based on the input NetCDF file. The ts_time variable
%           is just the input file times in Modified Julian Day.
%
% EXAMPLE USAGE
%    Mobj = get_AMM_forcing(Mobj, ts, depth)
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2013-02-07 First version based on get_POLCOMS_tsobc.m.
%
%==========================================================================

subname = 'get_AMM_tsobc';

global ftbverbose;
if ftbverbose
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end

varlist = {'lon', 'lat', 'ETWD', 'x1XD', 'time', 'depth', 'pdepthD'};

% Get the results. Check we have a cell array, and if we don't, assume it's
% a file name.

if iscell(ts)
    todo = length(ts);
else
    todo = 1;
end

for ii = 1:todo

    if ftbverbose
        if iscell(ts)
            ftn = ts{ii};
        else
            ftn = ts;
        end
        % Strip path from filename for the verbose output.
        [~, basename, ext] = fileparts(ftn);
        tmp_fn = [basename, ext];

        if todo == 1
            fprintf('%s: extracting file %s... ', subname, tmp_fn)
        else
            fprintf('%s: extracting file %s (%i of %i)... ', subname, tmp_fn, ii, todo)
        end
    end

    nc = netcdf.open(ftn, 'NOWRITE');

    for var = 1:numel(varlist)

        getVar = varlist{var};
        varid_amm = netcdf.inqVarID(nc, getVar);

        data = double(netcdf.getVar(nc, varid_amm, 'single'));
        if ii == 1
            amm.(getVar).data = data;
        else
            if ndims(data) < 4
                if strcmpi(varlist{var}, 'time')
                    % If the dimension is time, we need to be a bit more
                    % clever since we'll need a concatenated time series
                    % (in which values are continuous and from which we
                    % can extract a sensible time). As such, we need to add
                    % the maximum of the existing time. On the first
                    % iteration, we should save ourselves the base time
                    % (from the units of the time variable).
                    amm.(getVar).data = [amm.(getVar).data; data + max(amm.(getVar).data)];
                else
                    % This should be a fixed set of values (probably lon or
                    % lat) in which case we don't need to append them, so
                    % just replace the existing values with those in the
                    % current NetCDF file.
                    amm.(getVar).data = data;
                end
            elseif ndims(data) == 4
                % Assume we're concatenating with time (so along the fourth
                % dimesnion.
                amm.(getVar).data = cat(4, amm.(getVar).data, data);
            else
                error('Unsupported number of dimensions in AMM data')
            end
        end
        % Try to get some units (important for the calculation of MJD).
        try
            if ii == 1
                units = netcdf.getAtt(nc, varid_amm, 'units');
            else
                % Leave the units values alone so we always use the values
                % from the first file. This is particularly important for
                % the time calculation later on which is dependent on
                % knowing the time origin of the first file.
                continue
            end
        catch
            units = [];
        end
        amm.(getVar).units = units;
    end

    netcdf.close(nc)

    if ftbverbose
        fprintf('done.\n')
    end

end

% Data format:
% 
%   amm.ETWD.data and amm.x1XD.data are y, x, sigma, time
% 
[~, ~, nz, nt] = size(amm.ETWD.data);

% Make rectangular arrays for the nearest point lookup.
[lon, lat] = meshgrid(amm.lon.data, amm.lat.data);

fvlon = Mobj.lon(Mobj.obc_nodes(Mobj.obc_nodes ~= 0));
fvlat = Mobj.lat(Mobj.obc_nodes(Mobj.obc_nodes ~= 0));

% Number of boundary nodes
nf = sum(Mobj.nObcNodes);
% Number of sigma layers.
fz = size(Mobj.siglayz, 2);

fvtemp = nan(nf, fz, nt); % FVCOM interpolated temperatures
fvsal = nan(nf, fz, nt); % FVCOM interpolated salinities

if ftbverbose
    tic
end
for t = 1:nt
    % Get the current 3D array of AMM results.
    ammtemp3 = amm.ETWD.data(:, :, :, t);
    ammsalt3 = amm.x1XD.data(:, :, :, t);
    
    % Preallocate the intermediate results arrays.
    itempz = nan(nf, nz);
    isalz = nan(nf, nz);
    idepthz = nan(nf, nz);
    
    for j = 1:nz
        % Now extract the relevant layer from the 3D subsets. Transpose the
        % data to be (x, y) rather than (y, x).
        ammtemp2 = ammtemp3(:, :, j)';
        ammsalt2 = ammsalt3(:, :, j)';
        ammdepth2 = squeeze(amm.depth.data(:, :, j, t))';
       
        % Create new arrays which will be flattened when masking (below).
        tammtemp2 = ammtemp2;
        tammsalt2 = ammsalt2;
        tammdepth2 = ammdepth2;
        tlon = lon;
        tlat = lat;
        
        % Create and apply a mask to remove values outside the domain. This
        % inevitably flattens the arrays, but it shouldn't be a problem
        % since we'll be searching for the closest values in such a manner
        % as is appropriate for an unstructured grid (i.e. we're assuming
        % the AMM data is irregularly spaced).
        mask = tammdepth2 > 20000;
        tammtemp2(mask) = [];
        tammsalt2(mask) = [];
        tammdepth2(mask) = [];
        % Also apply the masks to the position arrays so we can't even find
        % positions outside the domain, effectively meaning if a value is
        % outside the domain, the nearest value to the boundary node will
        % be used.
        tlon(mask) = [];
        tlat(mask) = [];
        
        % Preallocate the intermediate results arrays.
        itempobc = nan(nf, 1);
        isalobc = nan(nf, 1);
        idepthobc = nan(nf, 1);
        
        % Speed up the tightest loop with a parallelized loop.
        parfor i = 1:nf
            % Now we can do each position within the 2D layer.

            fx = fvlon(i);
            fy = fvlat(i);

            [~, ii] = sort(sqrt((tlon - fx).^2 + (tlat - fy).^2));
            % Get the n nearest nodes from AMM (more? fewer?).
            ixy = ii(1:16);

            % Get the variables into static variables for the
            % parallelisation.
            plon = tlon(ixy);
            plat = tlat(ixy);
            ptemp = tammtemp2(ixy);
            psal = tammsalt2(ixy);
            pdepth = tammdepth2(ixy);
            
            % Use a triangulation to do the horizontal interpolation.
            tritemp = TriScatteredInterp(plon', plat', ptemp', 'natural');
            trisal = TriScatteredInterp(plon', plat', psal', 'natural');
            triz = TriScatteredInterp(plon', plat', pdepth', 'natural');
            itempobc(i) = tritemp(fx, fy);
            isalobc(i) = trisal(fx, fy);
            idepthobc(i) = triz(fx, fy);
            
            % Check all three, though if one is NaN, they all will be.
            if isnan(itempobc(i)) || isnan(isalobc(i)) || isnan(idepthobc(i))
                warning('FVCOM boundary node at %f, %f is outside the AMM domain. Setting to the closest AMM value.', fx, fy)
                itempobc(i) = tammtemp2(ii(1));
                isalobc(i) = tammsalt2(ii(1));
                idepthobc(i) = tammdepth2(ii(1));
            end
        end
        
        % Put the results in this intermediate array.
        itempz(:, j) = itempobc;
        isalz(:, j) = isalobc;
        idepthz(:, j) = idepthobc;
    end

    % Now we've interpolated in space, we can interpolate the z-values
    % to the sigma depths.
    oNodes = Mobj.obc_nodes(Mobj.obc_nodes ~= 0);
    for zi = 1:fz

        % Preallocate the output arrays
        fvtempz = nan(nf, fz);
        fvsalz = nan(nf, fz);

        for pp = 1:nf
            % Get the FVCOM depths at this node
            tfz = Mobj.siglayz(oNodes(pp), :);
            % Now get the interpolated AMM depth at this node
            tpz = idepthz(pp, :);

            % Get the temperature and salinity values for this node and
            % interpolate down the water column (from AMM to FVCOM). I had
            % originally planned to use csaps for the vertical
            % interplation/subsampling at each location. However, the demo
            % of csaps in the MATLAB documentation makes the interpolation
            % look horrible (shaving off extremes). interp1 provides a
            % range of interpolation schemes, of which pchip seems to do a
            % decent job of the interpolation (at least qualitatively).
            if ~isnan(tpz)
                fvtempz(pp, :) = interp1(tpz, itempz(pp, :), tfz, 'linear', 'extrap');
                fvsalz(pp, :) = interp1(tpz, isalz(pp, :), tfz, 'linear', 'extrap');
            else
                warning('Should never see this... ') % because we test for NaNs when fetching the values.
                warning('FVCOM boundary node at %f, %f is outside the AMM domain. Skipping.', fvlon(pp), fvlat(pp))
                continue
            end
        end
    end
    
    % The horizontally- and vertically-interpolated values in the final
    % FVCOM results array.
    fvtemp(:, :, t) = fvtempz;
    fvsal(:, :, t) = fvsalz;
end
if ftbverbose
    toc
end

Mobj.temperature = fvtemp;
Mobj.salt = fvsal;

% Convert the current times to Modified Julian Day (this is a bit ugly).
amm.time.all = strtrim(regexp(amm.time.units, 'since', 'split'));
amm.time.datetime = strtrim(regexp(amm.time.all{end}, ' ', 'split'));
amm.time.ymd = str2double(strtrim(regexp(amm.time.datetime{1}, '-', 'split')));
amm.time.hms = str2double(strtrim(regexp(amm.time.datetime{2}, ':', 'split')));

Mobj.ts_times = greg2mjulian(...
    amm.time.ymd(1), ...
    amm.time.ymd(2), ...
    amm.time.ymd(3), ...
    amm.time.hms(1), ...
    amm.time.hms(2), ...
    amm.time.hms(3)) + (amm.time.data / 3600 / 24);

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end