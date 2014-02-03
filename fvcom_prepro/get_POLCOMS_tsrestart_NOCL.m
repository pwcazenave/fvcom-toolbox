function restart = get_POLCOMS_tsrestart_NOCL(Mobj, inputConf)
% Extract temperature and salinity forcing information from NOC Operational
% Tide Surge model output and interpolate onto the FVCOM domain for use in 
% a reproduced FVCOM restart file.
%
% function restart = get_POLCOMS_tsrestart_NOCL(Mobj, inputConf)
%
% DESCRIPTION:
%    Interpolate temperature and salinity values onto the FVCOM nodes
%    at all sigma levels.
%
% INPUT:
%   Mobj        = MATLAB mesh structure which must contain:
%                   - Mobj.siglayz - sigma layer depths for all model
%                   nodes.
%                   - Mobj.lon, Mobj.lat and/or Mobj.x, Mobj.y - node
%                   coordinates.
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
% 
% OUTPUT:
%    restart = MATLAB structure which contains three fields (called
%              temperature, salinity and ts_time). temperature and salinity
%              have sizes (Mobj.nVerts, sigma, time). The time dimension is
%              determined based on the input NetCDF file. The ts_time
%              variable is just the input file times in Modified Julian Day.
%
% EXAMPLE USAGE
%    restart = get_POLCOMS_tsrestart_NOCL(Mobj, inputConf)
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Karen Amoudry (National Oceanography Centre, Liverpool)
%
% PWC Revision history
%    2013-01-09 First version based on the FVCOM shelf model
%    get_POLCOMS_forcing.m script (i.e. not a function but a plain script).
%
% KJA Revision history:
%    2014-01-15 First version, adapted from KJA's
%    'get_POLCOMS_tsobc_NOCL.m', which in turn was based on PWC's
%    'get_POLCOMS_tsobc.m'.
%
%==========================================================================

subname = 'get_POLCOMS_tsrestart_NOCL';

global ftbverbose;
if ftbverbose
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end
%%
% Which variables do we want from the POLCOMS file?
varlist = {'lon', 'lat', 'tem', 'sal', 'time'};

% Create the POLCOMS filename based on the year and month of interest
polcoms_ts = [inputConf.polcoms_ts,num2str(inputConf.startDate(1)),...
    '-',num2str(inputConf.startDate(2),'%02i'),'.nc'];

% Get the results from the POLCOMS file
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
% Convert FVCOM timestep into AMM/S12 output (seconds since 20071101:000000)
AMM_time = etime(inputConf.startDate,[2007,11,01,0,0,0]);
keep_time = (pc.time.data == AMM_time+(12*60*60)); % Allowance for POLCOMS data being at noon, not midnight
pc.tem.data = pc.tem.data(:,:,:,keep_time);
pc.sal.data = pc.sal.data(:,:,:,keep_time);

% Import the POLCOMS sigma info
pc = get_POLCOMS_sigma(pc,inputConf);

% Data format:
%   pc.tem.data and pc.sal.data are x, y, sigma
[~, ~, nz] = size(pc.tem.data);

% Make rectangular arrays for the nearest point lookup.
[plon, plat] = meshgrid(pc.lon.data, pc.lat.data);

% Number of sigma layers.
fz = size(Mobj.siglayz, 2);

if ftbverbose
    tic
end

%  Flip the vertical layer dimension to make the POLCOMS data go from
% surface to seabed to match its depth data and to match how FVCOM
% works.
pctemp3 = flipdim(pc.tem.data, 3);
pcsalt3 = flipdim(pc.sal.data, 3);

% Preallocate the FVCOM results arrays
fvtemp = nan(Mobj.nVerts, fz); % FVCOM interpolated temperatures
fvsal = nan(Mobj.nVerts, fz); % FVCOM interpolated salinities

% Preallocate the intermediate results arrays.
itempz = nan(Mobj.nVerts, nz);
isalz = nan(Mobj.nVerts, nz);
idepthz = nan(Mobj.nVerts,nz);

for j = 1:nz
    % Transpose the POLCOMS data to be (x,y) rather than (y,x)
    pctemp2 = pctemp3(:, :, j)';
    pcsalt2 = pcsalt3(:, :, j)';
    pcdepth2 = squeeze(pc.depth.data(:, :, j))';
    
    % Reshape the arrays to allow the sort to work properly later
    tlon=reshape(plon,(size(plon,1)*size(plon,2)),1);
    tlat=reshape(plat,(size(plat,1)*size(plat,2)),1);
    pctemp2 = reshape(pctemp2,(size(pctemp2,1)*size(pctemp2,2)),1);
    pcsalt2 = reshape(pcsalt2,(size(pcsalt2,1)*size(pcsalt2,2)),1);
    pcdepth2 = reshape(pcdepth2,(size(pcdepth2,1)*size(pcdepth2,2)),1);
    
    % Find the points which aren't NaNs
    keeptemp = find(~isnan(pctemp2));
    keepsalt = find(~isnan(pcsalt2));
    keepdepth = find(~isnan(pcdepth2));
    
    keep = intersect(keeptemp,keepsalt);
    keep = intersect(keep, keepdepth);
    
    % Interpolate the POLCOMS results by (1) turning the non-NaN POLCOMS
    % results into a Tri object and (2) extracting the values of that Tri
    % Object at the FVCOM node locations.
    tritemp = TriScatteredInterp(tlon(keep),tlat(keep),pctemp2(keep),'natural');
    itemp = tritemp(Mobj.lon,Mobj.lat);
    
    trisalt = TriScatteredInterp(tlon(keep),tlat(keep),pcsalt2(keep),'natural');
    isalt = trisalt(Mobj.lon,Mobj.lat);
    
    tridepth = TriScatteredInterp(tlon(keep),tlat(keep),pcdepth2(keep),'natural');
    idepth = tridepth(Mobj.lon,Mobj.lat);
    
    % Put the results in this intermediate array.
    itempz(:, j) = itemp;
    isalz(:, j) = isalt;
    idepthz(:,j) = idepth;
end

% Now we've interpolated in space, we can interpolate the z-values
% to the sigma depths.
% for zi = 1:fz
    for pp = 1:Mobj.nVerts
        % Get the FVCOM depths at this node
        tfz = Mobj.siglayz(pp, :);
        % Now get the interpolated POLCOMS depth at this node
        tpz = idepthz(pp, :);
        
        % To ensure we get the full vertical expression of the vertical
        % profiles, we need to normalise the POLCOMS and FVCOM
        % depths to the same range. This is because in instances where
        % FVCOM depths are shallower (e.g. in coastal regions), if we
        % don't normalise the depths, we end up truncating the vertical
        % profile. This approach ensures we always use the full
        % vertical profile, but we're potentially squeezing it into a
        % smaller depth.
        A = max(tpz);
        B = min(tpz);
        C = max(tfz);
        D = min(tfz);
        norm_tpz = (((D - C) * (tpz - A)) / (B - A)) + C;
        
        % Get the temperature and salinity values for this node and
        % interpolate down the water column (from POLCOMS to FVCOM).
        % Change to 'pchip' to match PWC parent code.
        if ~isnan(tpz)
            fvtemp(pp, :) = interp1(norm_tpz, itempz(pp, :), tfz, 'pchip', 'extrap');
            fvsal(pp, :) = interp1(norm_tpz, isalz(pp, :), tfz, 'pchip', 'extrap');
        else
            warning('Should never see this... ') % because we test for NaNs when fetching the values.
            warning('FVCOM boundary node at %f, %f is outside the POLCOMS domain. Skipping.', Mobj.lon(pp), Mobj.lat(pp))
            continue
        end
    end
% end

if ftbverbose
    toc
end

% Convert NOC model output temperatures from Kelvin to Celsius
fvtemp = fvtemp - 273.15;

% Timeshift to match the expected FVCOM input times. The temperature and
% salinity values are a daily average (midnight to midnight). They are
% given a timestamp of the middle time in this period (noon). Each day's
% value applies to the whole day (e.g. from midnight to midnight).
% Therefore, we can shift the time to midnight at the start of the relevant
% day and apply the value to the whole day. This makes FVCOM (and me)
% happy.
restart.ts_times = datenum(2007, 11, 0, 0, 0, 0) + ...    % Convert NOC model reference time to Matlab datenum
    ((pc.time.data(keep_time)+(12*3600)) / 3600 / 24);  % Add NOC model output time and 12 hour offset
restart.ts_times = datevec(restart.ts_times);   % Convert times to datevec format

% Final output variable names must match FVCOM restart file variable names
restart.temp = fvtemp;
restart.salinity = fvsal;

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end