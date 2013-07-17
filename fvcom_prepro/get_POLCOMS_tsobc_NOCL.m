function Mobj = get_POLCOMS_tsobc_NOCL(Mobj, inputConf)
% Extract temperature and salinity boundary forcing information from NOC
% Operation Tide Surge model output.
%
% function Mobj = get_POLCOMS_tsobc_NOCL(Mobj, inputConf)
%
% DESCRIPTION:
%    Interpolate temperature and salinity values onto the FVCOM open
%    boundaries at all sigma levels.
%
% INPUT:
%   Mobj        = MATLAB mesh structure which must contain:
%                   - Mobj.siglayz - sigma layer depths for all model
%                   nodes.
%                   - Mobj.lon, Mobj.lat and/or Mobj.x, Mobj.y - node
%                   coordinates.
%                   - Mobj.obc_nodes - list of open boundary node inidices.
%                   - Mobj.nObcNodes - number of nodes in each open
%                   boundary.
%   inputConf   = MATLAB structure which must contain: 
%                   - inputConf.polcoms_ts - location of NOC Operational
%                   Model output containing 4D variables of temperature
%                   (tem) and salinity (sal). They should have dimensions
%                   (x, y, sigma, time).
%                   - inputConf.polcoms_z - location of NOC Operational
%                   Model output containing 4D variables of bathymetry
%                   (XXX) and sigma layer thickness (XXX).
%                   - inputConf.startDate - start date and time for FVCOM
%                   model run
%                   - inputConf.endDate - end date and time for FVCOM
%                   model run
% 
% OUTPUT:
%    Mobj = MATLAB structure in which three new fields (called temperature,
%           salinity and ts_time). temperature and salinity have sizes
%           (sum(Mobj.nObcNodes), sigma, time). The time dimension is
%           determined based on the input NetCDF file. The ts_time variable
%           is just the input file times in Modified Julian Day.
%
% EXAMPLE USAGE
%    Mobj = get_POLCOMS_forcing_NOCL(Mobj, inputConf)
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Karen Thurston (National Oceanography Centre, Liverpool)
%
% PWC Revision history
%    2013-01-09 First version based on the FVCOM shelf model
%    get_POLCOMS_forcing.m script (i.e. not a function but a plain script).
%
% KJT Revision history:
%    2013-02-05 Adapted from PWC's script to fit NOCL file formats.
%
%==========================================================================

subname = 'get_POLCOMS_tsobc_NOCL';

global ftbverbose;
if ftbverbose
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end
%%
varlist = {'lon', 'lat', 'tem', 'sal', 'time'};

if (inputConf.startDate(:,2)==inputConf.startDate(:,2))
    polcoms_ts = [inputConf.polcoms_ts,num2str(inputConf.startDate(:,1)),...
        '-',num2str(inputConf.startDate(:,2),'%02i'),'.nc'];
else
    error('Model run spans two months. I have not prepared for this eventuality yet.')
end

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

% Get rid of data outside the time we're interested in
% Create an array of hourly timesteps, ensuring the output time series is
% at least as long as the FVCOM model run time.
timesteps = datevec(datenum(inputConf.startDate):1/24:datenum(inputConf.endDate)+1);

% Convert FVCOM timestep into AMM/S12 output (seconds since 20071101:000000)
AMM_time = etime(timesteps,repmat([2007,11,01,0,0,0],size(timesteps,1),1));
keep_time = (pc.time.data >= AMM_time(1) & pc.time.data <= AMM_time(end));
pc.tem.data = pc.tem.data(:,:,:,keep_time);
pc.sal.data = pc.sal.data(:,:,:,keep_time);

pc = get_POLCOMS_sigma(pc,inputConf);

% Data format:
% 
%   pc.tem.data and pc.sal.data are x, y, sigma, time
% 
[~, ~, nz, nt] = size(pc.tem.data);

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
    pctemp2 = pc.tem.data(:, :, :, t);
    pcsal2 = pc.sal.data(:, :, :, t);
    
    % Preallocate the intermediate results arrays.
    itempz = nan(nf, nz);
    isalz = nan(nf, nz);
    idepthz = nan(nf,nz);
    
    for j = 1:nz
        pcdepth2 = squeeze(pc.depth.data(:, :, j));
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
        % Reshape the arrays to allow the sort to work properly later
        tlon=reshape(tlon,(size(tlon,1)*size(tlon,2)),1);
        tlat=reshape(tlat,(size(tlat,1)*size(tlat,2)),1);
        
        % Preallocate the intermediate results arrays.
        itempobc = nan(nf, 1);
        isalobc = nan(nf, 1);
        idepthobc = nan(nf, 1);
        
        % Speed up the tightest loop with a parallelized loop.
        for i = 1:nf
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
            tritemp = TriScatteredInterp(plon, plat, ptemp, 'natural');
            trisal = TriScatteredInterp(plon, plat, psal, 'natural');
            triz = TriScatteredInterp(plon, plat, pdepth, 'natural');
            itempobc(i) = tritemp(fx, fy);
            isalobc(i) = trisal(fx, fy);
            idepthobc(i) = triz(fx, fy);
            
            % Check all three, though if one is NaN, they all will be.
            if isnan(itempobc(i)) || isnan(isalobc(i)) || isnan(idepthobc(i))
                if ftbverbose
                    warning('FVCOM boundary node at %f, %f is outside the POLCOMS domain. Setting to the closest POLCOMS value.', fx, fy)
                end
                p = 1;
                while isnan(tpctemp2(ii(p)))
                    p = p+1;
                end
                itempobc(i) = tpctemp2(ii(p));
                p = 1;
                while isnan(tpcsal2(ii(p)))
                    p = p+1;
                end
                isalobc(i) = tpcsal2(ii(p));
                p = 1;
                while isnan(tpcdepth2(ii(p)))
                    p = p+1;
                end                          
                idepthobc(i) = tpcdepth2(ii(p));
            end
        end
        
        % Put the results in this intermediate array.
        itempz(:, j) = itempobc;
        isalz(:, j) = isalobc;
        idepthz(:,j) = idepthobc;
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
            % each location.
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

% Convert NOC model output temperatures from Kelvin to Celsius
fvtemp = fvtemp - 273.15;

Mobj.temperature = fvtemp;
Mobj.salt = fvsal;

% Interpolate to the FVCOM time series. Operational model daily output is
% at noon, and we want midnight.
Mobj.ts_times = greg2mjulian(2007, 11, 0, 0, 0, 0) + ...    % Convert NOC model reference time to MJD
    ((pc.time.data(keep_time)+(12*3600)) / 3600 / 24);  % Add NOC model output time and 12 hour offset
 

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end