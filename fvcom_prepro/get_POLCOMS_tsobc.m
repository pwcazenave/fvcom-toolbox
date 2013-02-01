function Mobj = get_POLCOMS_tsobc(Mobj, polcoms_ts, polcoms_z)
% Read temperature and salinity from POLCOMS NetCDF model output files and
% interpolate onto the open boundaries in Mobj.
%
% function Mobj = get_POLCOMS_tsobc(Mobj, polcoms_ts, polcoms_bathy, varlist)
%
% DESCRIPTION:
%    Interpolate temperature and salinity values onto the FVCOM open
%    boundaries at all sigma levels.
%
% INPUT:
%   Mobj        = MATLAB mesh structure which must contain:
%                   - Mobj.siglayz - sigma layer depths for all model
%                   nodes.
%                   - Mobj.lon, Mobj.lat - node coordinates (lat/long).
%                   - Mobj.obc_nodes - list of open boundary node inidices.
%                   - Mobj.nObcNodes - number of nodes in each open
%                   boundary.
%   polcoms_ts  = POLCOMS NetCDF file in which 4D variables of temperature 
%                 and salinity (called 'ETW' and 'x1X') exist. Their shape
%                 should be (y, x, sigma, time).
%   polcoms_z   = POLCOMS NetCDF file in which 4D variables of bathymetry
%                 and sigma layer thickness can be found. They should be
%                 called 'depth' and 'pdepth' respectively.
% 
% OUTPUT:
%    Mobj = MATLAB structure in which three new fields (called temperature,
%           salinity and ts_time). temperature and salinity have sizes
%           (sum(Mobj.nObcNodes), sigma, time). The time dimension is
%           determined based on the input NetCDF file. The ts_time variable
%           is just the input file times in Modified Julian Day.
%
% EXAMPLE USAGE
%    Mobj = get_POLCOMS_forcing(Mobj, polcoms_ts, polcoms_z)
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2013-01-09 First version based on the FVCOM shelf model
%    get_POLCOMS_forcing.m script (i.e. not a function but a plain script).
%
%==========================================================================

subname = 'get_POLCOMS_forcing';

global ftbverbose;
if ftbverbose
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end

varlist = {'lon', 'lat', 'ETW', 'x1X', 'time'};

% Get the results
nc = netcdf.open(polcoms_ts, 'NOWRITE');

for var=1:numel(varlist)
    
    getVar = varlist{var};
    varid_pc = netcdf.inqVarID(nc, getVar);
    
    data = netcdf.getVar(nc, varid_pc, 'single');
    pc.(getVar).data = double(data);
    % Try to get some units (important for the calculation of MJD).
    try
        units = netcdf.getAtt(nc,varid_pc,'units');
    catch
        units = [];
    end
    pc.(getVar).units = units;
end

netcdf.close(nc)

% Extract the bathymetry ('pdepth' is cell thickness, 'depth' is cell
% centre depth).
nc = netcdf.open(polcoms_z, 'NOWRITE');
varid_pc = netcdf.inqVarID(nc, 'depth'); 
pc.depth.data = double(netcdf.getVar(nc, varid_pc, 'single'));
try
    pc.depth.units = netcdf.getAtt(nc, varid_pc, 'units');
catch
    pc.depth.units = [];
end
varid_pc = netcdf.inqVarID(nc, 'pdepth'); 
pc.pdepth.data = double(netcdf.getVar(nc, varid_pc, 'single'));
try
    pc.pdepth.units = netcdf.getAtt(nc, varid_pc, 'units');
catch
    pc.pdepth.units = [];
end
netcdf.close(nc)

% Data format:
% 
%   pc.ETW.data and pc.x1X.data are y, x, sigma, time
% 
[~, ~, nz, nt] = size(pc.ETW.data);

% Make rectangular arrays for the nearest point lookup.
[lon, lat] = meshgrid(pc.lon.data, pc.lat.data);

% Find the nearest POLCOMS point to each point in the FVCOM open boundaries
fvlon = Mobj.lon(Mobj.obc_nodes(Mobj.obc_nodes ~= 0));
fvlat = Mobj.lat(Mobj.obc_nodes(Mobj.obc_nodes ~= 0));

% Number of boundary nodes
nf = sum(Mobj.nObcNodes);
% Number of sigma layers.
fz = size(Mobj.siglayz, 2);

% itemp = nan(nf, nz, nt); % POLCOMS interpolated temperatures
% isal = nan(nf, nz, nt); % POLCOMS interpolated salinities
fvtemp = nan(nf, fz, nt); % FVCOM interpolated temperatures
fvsal = nan(nf, fz, nt); % FVCOM interpolated salinities

if ftbverbose
    tic
end
for t = 1:nt
    % Get the current 3D array of POLCOMS results.
    pctemp3 = pc.ETW.data(:, :, :, t);
    pcsal3 = pc.x1X.data(:, :, :, t);
    
    % Preallocate the intermediate results arrays.
    itempz = nan(nf, nz);
    isalz = nan(nf, nz);
    idepthz = nan(nf, nz);
    
    for j = 1:nz
        % Now extract the relevant layer from the 3D subsets. Transpose the
        % data to be (x, y) rather than (y, x).
        pctemp2 = pctemp3(:, :, j)';
        pcsal2 = pcsal3(:, :, j)';
        pcdepth2 = squeeze(pc.depth.data(:, :, j, t))';
       
        % Create new arrays which will be flattened when masking (below).
        tpctemp2 = pctemp2;
        tpcsal2 = pcsal2;
        tpcdepth2 = pcdepth2;
        tlon = lon;
        tlat = lat;
        
        % Create and apply a mask to remove values outside the domain. This
        % inevitably flattens the arrays, but it shouldn't be a problem
        % since we'll be searching for the closest values in such a manner
        % as is appropriate for an unstructured grid (i.e. we're assuming
        % the POLCOMS data is irregularly spaced).
        mask = tpcdepth2 < -20000;
        tpctemp2(mask) = [];
        tpcsal2(mask) = [];
        tpcdepth2(mask) = [];
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
            % Get the n nearest nodes from POLCOMS (more? fewer?).
            ixy = ii(1:16);

            % Get the variables into static variables for the
            % parallelisation.
            plon = tlon(ixy);
            plat = tlat(ixy);
            ptemp = tpctemp2(ixy);
            psal = tpcsal2(ixy);
            pdepth = tpcdepth2(ixy);
            
            % Use a triangulation to do the horizontal interpolation.
            tritemp = TriScatteredInterp(plon', plat', ptemp', 'natural');
            trisal = TriScatteredInterp(plon', plat', psal', 'natural');
            triz = TriScatteredInterp(plon', plat', pdepth', 'natural');
            itempobc(i) = tritemp(fx, fy);
            isalobc(i) = trisal(fx, fy);
            idepthobc(i) = triz(fx, fy);
            
            % Check all three, though if one is NaN, they all will be.
            if isnan(itempobc(i)) || isnan(isalobc(i)) || isnan(idepthobc(i))
                warning('FVCOM boundary node at %f, %f is outside the POLCOMS domain. Setting to the closest POLCOMS value.', fx, fy)
                itempobc(i) = tpctemp2(ii(1));
                isalobc(i) = tpcsal2(ii(1));
                idepthobc(i) = tpcdepth2(ii(1));
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
            % Now get the interpolated POLCOMS depth at this node
            tpz = idepthz(pp, :);

            % Get the temperature and salinity values for this node and
            % interpolate down the water column (from POLCOMS to FVCOM).
            % TODO: Use csaps for the vertical interplation/subsampling at
            % each location. Alternatively, the pchip interp1 method seems
            % to do a decent job of the interpolation; it might be a more
            % suitable candidate in the absence of csaps. In fact, the demo
            % of csaps in the MATLAB documentation makes the interpolation
            % look horrible (shaving off extremes). I think pchip is
            % better.
            if ~isnan(tpz)
                fvtempz(pp, :) = interp1(tpz, itempz(pp, :), tfz, 'linear', 'extrap');
                fvsalz(pp, :) = interp1(tpz, isalz(pp, :), tfz, 'linear', 'extrap');
            else
                warning('Should never see this... ') % because we test for NaNs when fetching the values.
                warning('FVCOM boundary node at %f, %f is outside the POLCOMS domain. Skipping.', fvlon(pp), fvlat(pp))
                continue
            end
        end
    end
    
    % The horizontally-interpolated values in the final results array.
%     itemp(:, :, t) = itempz;
%     isal(:, :, t) = isalz;
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

% Do we have to interpolate to the FVCOM time series? Looking at page 325
% of the FVCOM manual, it looks like the temperature and salinity are on a
% different time sampling from the other example files (14 time steps vs.
% 3625 for _wnd.nc or 43922 for _julain_obc.nc (i.e. surface elevation at
% the boundary)). That's not to say those files were all used in the same
% model run... In the interim, just convert the current times to Modified
% Julian Day (this is a bit ugly).
% pc.time.yyyymmdd = strtrim(regexp(pc.time.units, 'since', 'split'));
% pc.time.yyyymmdd = str2double(regexp(pc.time.yyyymmdd{end}, '-', 'split'));
% Mobj.ts_times = greg2mjulian(pc.time.yyyymmdd(1), pc.time.yyyymmdd(2), pc.time.yyyymmdd(3), 0, 0, 0) + (pc.time.data / 3600 / 24);

% Convert the POLCOMS times to Modified Julian Day (this is a very ugly).
pc.time.yyyymmdd = strtrim(regexp(pc.time.units, 'since', 'split'));
pc.time.strtime = regexp(pc.time.yyyymmdd{end}, '-', 'split');
% This new version of the time has the year in a weird format (yr.#). We
% thus need to split it again to get the decimal year (post-2000 only?).
pc.time.strtimeyr = regexp(pc.time.strtime, '\.', 'split');
pc.time.yyyymmdd = str2double([pc.time.strtimeyr{1}(2), pc.time.strtime{2:3}]);
pc.time.yyyymmdd(1) = pc.time.yyyymmdd(1) + 2000; % add full year.
Mobj.ts_times = greg2mjulian(pc.time.yyyymmdd(1), pc.time.yyyymmdd(2), pc.time.yyyymmdd(3), 0, 0, 0) + (pc.time.data / 3600 / 24);

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end