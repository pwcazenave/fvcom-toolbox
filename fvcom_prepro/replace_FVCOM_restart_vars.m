function replace_FVCOM_restart_vars(Mobj, polcoms_ts, polcoms_z, start_date, fv_restart, varlist)
% Use an FVCOM restart file to seed a model run with spatially varying
% versions of otherwise constant variables (temperature and salinity only
% for the time being).
%
% function replace_FVCOM_restart_vars(Mobj, polcoms_ts, polcoms_z, ...
%    start_date, fv_restart, varlist)
%
% DESCRIPTION:
%    FVCOM does not yet support spatially varying temperature and salinity
%    inputs as initial conditions. To avoid having to run a model for a
%    long time in order for temperature and salinity to settle within the
%    model from the atmospheric and boundary forcing, we can use a restart
%    file to cheat. If we run a model for a week or so (until the
%    hydrodynamics has stabilised, we can use the restart file generated
%    from that run as the basis for a new run, except we replace the
%    currently computed temperature and salinity and replace them with data
%    interpolated from another source (in this case, POLCOMS). 
%
% INPUT:
%   Mobj        = MATLAB mesh structure which must contain:
%                   - Mobj.siglayz - sigma layer depths for all model
%                   nodes.
%                   - Mobj.lon, Mobj.lat - node coordinates (long/lat).
%   polcoms_ts  = POLCOMS NetCDF file in which 4D variables of temperature
%   and salinity (called 'ETW' and 'x1X') exist. Their shape should be (y,
%   x, sigma, time).
%   polcoms_z   = POLCOMS NetCDF file in which 4D variables of bathymetry
%   and sigma layer thickness can be found. They should be called 'depth'
%   and 'pdepth' respectively.
%   start_date  = Gregorian start date array (YYYY, MM, DD, hh, mm, ss).
%   fv_restart  = full path to the FVCOM restart file.
%   varlist     = cell array of variables to extract from the POLCOMS
%   NetCDF files (polcoms_ts only, pdepth and depth will be extracted from
%   polcoms_z).
% 
% OUTPUT:
%   NetCDF file in which the temperature and salinity variables have been
%   replaced with the POLCOMS versions. The file name is the input restart
%   file name appended _polcoms.nc.
%
% EXAMPLE USAGE
%   replace_FVCOM_restart_vars(Mobj, '/tmp/pc_ts.nc', '/tmp/pc_z.nc', ...
%       '2006-01-01 00:00:00', '/tmp/fvcom_restart.nc', ...
%       {'lon', 'lat', 'ETW', 'x1X', 'time'})
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2013-01-28 First version (based on initial_tsobc.m).
%
%==========================================================================

subname = 'replace_FVCOM_restart_vars';

global ftbverbose;
if ftbverbose
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end

%--------------------------------------------------------------------------
% Extract the POLCOMS data specified in varlist
%--------------------------------------------------------------------------

if ftbverbose
    fprintf('%s : read POLCOMS data... ', subname)
end

nc = netcdf.open(polcoms_ts, 'NOWRITE');

for var=1:numel(varlist)
    
    getVar = varlist{var};
    varid_pc = netcdf.inqVarID(nc, getVar);
    
    data = netcdf.getVar(nc, varid_pc, 'single');
    pc.(getVar).data = double(data);
    % Try to get some units (important for the calculation of MJD).
    try
        units = netcdf.getAtt(nc,varid_pc,'units');
    catch err
        units = [];
    end
    pc.(getVar).units = units;
end
clear data getVar varid_pc var

netcdf.close(nc)

% Extract the bathymetry ('pdepth' is cell thickness, 'depth' is cell
% centre depth).
nc = netcdf.open(polcoms_z, 'NOWRITE');
varid_pc = netcdf.inqVarID(nc, 'depth'); 
pc.depth.data = double(netcdf.getVar(nc, varid_pc, 'single'));
try
    pc.depth.units = netcdf.getAtt(nc, varid_pc, 'units');
catch err
    pc.depth.units = [];
end
varid_pc = netcdf.inqVarID(nc, 'pdepth'); 
pc.pdepth.data = double(netcdf.getVar(nc, varid_pc, 'single'));
try
    pc.pdepth.units = netcdf.getAtt(nc, varid_pc, 'units');
catch err
    pc.pdepth.units = [];
end
netcdf.close(nc)

% Data format:
% 
%   pc.ETW.data and pc.x1X.data are y, x, sigma, time
% 
[ny, nx, ~, ~] = size(pc.ETW.data);

% Number of sigma layers.
[fn, fz] = size(Mobj.siglayz);

% Make rectangular arrays for the nearest point lookup.
[lon, lat] = meshgrid(pc.lon.data, pc.lat.data);

% Convert the POLCOMS times to Modified Julian Day (this is a very ugly).
pc.time.yyyymmdd = strtrim(regexp(pc.time.units, 'since', 'split'));
pc.time.strtime = regexp(pc.time.yyyymmdd{end}, '-', 'split');
% This new version of the time has the year in a weird format (yr.#). We
% thus need to split it again to get the decimal year (post-2000 only?).
pc.time.strtimeyr = regexp(pc.time.strtime, '\.', 'split');
pc.time.yyyymmdd = str2double([pc.time.strtimeyr{1}(2), pc.time.strtime{2:3}]);
pc.time.yyyymmdd(1) = pc.time.yyyymmdd(1) + 2000; % add full year.
Mobj.ts_times = greg2mjulian(pc.time.yyyymmdd(1), pc.time.yyyymmdd(2), pc.time.yyyymmdd(3), 0, 0, 0) + (pc.time.data / 3600 / 24);

% Given our intput time (in start_date), find the nearest time
% index for the POLCOMS data.
stime = greg2mjulian(start_date(1), start_date(2), ...
    start_date(3), start_date(4), ...
    start_date(5), start_date(6));
[~, tidx] = min(abs(Mobj.ts_times - stime));

if ftbverbose
    fprintf('done.\n') 
end

%--------------------------------------------------------------------------
% Interpolate the regularly gridded POLCOMS data onto the FVCOM grid
% (vertical grid first).
%--------------------------------------------------------------------------

if ftbverbose
    fprintf('%s : interpolate POLCOMS onto FVCOM''s vertical grid... ', subname)
end

% Preallocate the output arrays
pctempz = nan(nx, ny, fz);
pcsalz = nan(nx, ny, fz);

for xi = 1:nx
    % For the parallel loop, get all the y and z dimension data for the
    % current x position (temperature, salinity and depth).
    % N.B. The POLCOMS data is stored y, x, z, t (first two dimensions the
    % wrong way around).
    xtemp = squeeze(pc.ETW.data(:, xi, :, tidx));
    xsalt = squeeze(pc.x1X.data(:, xi, :, tidx));
    xdepth = squeeze(pc.depth.data(:, xi, :, tidx));
    
    % Preallocate the arrays for the inner loop
    ytemp = nan(ny, fz);
    ysalt = nan(ny, fz);
    for yi = 1:ny
        % First things first, check we're within the POLCOMS domain
        % (assuming anything beyond a maximum depth of 20km is outside the
        % domain).
        if xdepth(yi, end) < -20000 % use deepest depth
            continue
        end

        % Find the nearest sigma layer z values from the unstructured grid.
        [md, mi] = min(sqrt((Mobj.lon - lon(xi, yi)).^2 + (Mobj.lat - lat(xi, yi)).^2));

        % Skip data point if the closest FVCOM node is more than 10 minutes
        % away.
        if md > 10 / 60
            continue
        else
            % Use the FVCOM node's sigma depths to interpolate this POLCOMS
            % position's temperature and salinity data.
            
            % Get the FVCOM depths closest to this POLCOMS grid position.
            tfz = Mobj.siglayz(mi, :);
            % Now get the POLCOMS depth at this node for the time index
            % identified above.
            tpz = xdepth(yi, :); 
            ytemp(yi, :) = interp1(tpz, xtemp(yi, :), tfz, 'linear', 'extrap');
            
            ysalt(yi, :) = interp1(tpz, xsalt(yi, :), tfz, 'linear', 'extrap');
        end
    end
    pctempz(xi, :, :) = ytemp;
    pcsalz(xi, :, :) = ysalt;
end
% Tidy up the namespace a bit.
clear ytemp ysalt tfz tpz md mi xtemp xsalt xdepth xi yi zi

if ftbverbose
    fprintf('done.\n') 
end

%--------------------------------------------------------------------------
% Now we have vertically interpolated POLCOMS data, we can interpolate each
% sigma layer onto the FVCOM unstructured grid ready to write out to
% NetCDF. We'll use the triangular interpolation in MATLAB with the natural
% method (gives pretty good results, at least qualitatively).
%--------------------------------------------------------------------------

if ftbverbose
    fprintf('%s : interpolate POLCOMS onto FVCOM''s horizontal grid... ', subname)
end

fvtemp = nan(fn, fz);
fvsalt = nan(fn, fz);
for zi = 1:fz
    % Set up the interpolation object.
    ft = TriScatteredInterp(lon(:), lat(:), reshape(pctempz(:, :, zi), [], 1), 'natural');
    fs = TriScatteredInterp(lon(:), lat(:), reshape(pcsalz(:, :, zi), [], 1), 'natural');
    % Interpolate temperature and salinity onto the unstructured grid.
    fvtemp(:, zi) = ft(Mobj.lon, Mobj.lat);
    fvsalt(:, zi) = fs(Mobj.lon, Mobj.lat);
end

% Unfortunately, TriScatteredInterp won't extrapolate, returning instead
% NaNs outside the original data's extents. So, for each NaN position, find
% the nearest non-NaN value and use that instead. The order in which the
% NaN-nodes are found will determine the spatial pattern of the
% extrapolation.

% We can assume that all layers will have NaNs in the same place
% (horizontally), so just use the surface layer (1) for the identification
% of NaNs. Also store the finite values so we can find the nearest real
% value to the current NaN node and use its temperature and salinity
% values.
fvidx = 1:fn;
fvnanidx = fvidx(isnan(fvtemp(:, 1)));
fvfinidx = fvidx(~isnan(fvtemp(:, 1)));
for ni = 1:length(fvnanidx);
    % Current position
    xx = Mobj.lon(fvnanidx(ni));
    yy = Mobj.lat(fvnanidx(ni));
    % Find the nearest non-nan temperature and salinity value.
    [~, di] = min(sqrt((Mobj.lon(fvfinidx) - xx).^2 + (Mobj.lat(fvfinidx) - yy).^2));
    % Replace the temperature and salinity values at all depths at the
    % current NaN position with the closest non-nan value.
    fvtemp(fvnanidx(ni), :) = fvtemp(fvfinidx(di), :);
    fvsalt(fvnanidx(ni), :) = fvsalt(fvfinidx(di), :);
end

if ftbverbose
    fprintf('done.\n') 
end

%--------------------------------------------------------------------------
% Get the restart data and replace with the interpolated data.
%--------------------------------------------------------------------------

if ftbverbose
    fprintf('%s : export interpolated data to NetCDF:\n', subname)
end

nc = netcdf.open(fv_restart, 'NOWRITE');
[pp, nn, ee] = fileparts(fv_restart);
ncout = netcdf.create(fullfile(pp, [nn, '_polcoms', ee]), 'clobber');
% ncout = netcdf.create(fullfile('/tmp/pica/', [nn, '_polcoms', ee]), 'clobber');

[numdims, numvars, numglobatts, unlimdimID] = netcdf.inq(nc);

% Define the dimensions for all variables.
dimid = nan(numdims, 1);
dimnames = cell(numdims, 1);
dimlengths = nan(numdims, 1);
for ii = 1:numdims
    [dimname, dimlen] = netcdf.inqDim(nc, ii - 1);
    if ii ~= unlimdimID + 1 % NetCDF indices start at zero
        dimid(ii) = netcdf.defDim(ncout, dimname, dimlen);
    else
        dimid(ii) = netcdf.defDim(ncout, dimname, netcdf.getConstant('NC_UNLIMITED'));
    end
    dimnames{ii} = dimname;
    dimlengths(ii) = dimlen;
end

% Now define the variables and attributes.
for ii = 1:numvars

    % Find name of the current variable.
    [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(nc, ii - 1);

    % Create the variables.
    varid = netcdf.defVar(ncout, varname, xtype, varDimIDs);

    % Get each attribute and add it to the current variable.
    for j = 1:varAtts

        attname = netcdf.inqAttName(nc, varid, j - 1);
        attval = netcdf.getAtt(nc, varid, attname);

        netcdf.putAtt(ncout, varid, attname, attval);
    end
end

% Do the global attributes
for ii = 1:numglobatts
    
    % Find the current global attribute's name and value.
    gattname = netcdf.inqAttName(nc, netcdf.getConstant('NC_GLOBAL'), ii - 1);
    gattval = netcdf.getAtt(nc, netcdf.getConstant('NC_GLOBAL'), gattname);
    
    % Put that back into the output NetCDF file.
    netcdf.putAtt(ncout, netcdf.getConstant('NC_GLOBAL'), gattname, gattval);
end

netcdf.endDef(ncout);

% Get the existing data and output to the new NetCDF file, except for
% temperature and salinity, where we replace the values with the
% interpolated POLCOMS data.
for ii = 1:numvars

    [varname, ~, varDimIDs, ~] = netcdf.inqVar(nc, ii - 1);
    varid = netcdf.inqVarID(nc, varname);

    if ftbverbose
        fprintf('\tvariable %s... ', varname)
    end

    % We need the data irrespective of whether we're replacing it or not,
    % so grab it outside the if statement below.
    data = netcdf.getVar(nc, varid);

    % Get the size of the data and the dimension names.
    currDimsNames = dimnames(varDimIDs + 1);
    currDimsLengths = dimlengths(varDimIDs + 1);

    % Find whether we've got an unlimited dimension in this data.
    wasUnlimited = -1;
    for jj = varDimIDs
        if numel(unlimdimID) > 1
            error('Do not currently support multiple unlimited dimensions.')
        end
        if strcmpi(dimnames(jj + 1), dimnames(unlimdimID + 1))
            wasUnlimited = jj;
        end
    end

    % Since the restart file has a number of time values, we'll ramp up
    % the temperature from some constant to the actual value over the
    % time steps. So, we need to know how many time steps we actually
    % have.

    % Get the dimension data ready for the temperature and salinity arrays.
    tIdx = strncmp(dimnames(unlimdimID + 1), currDimsNames, length(dimnames{unlimdimID + 1}));
    % Not sure about the hardcoded strings below...
    sIdx = strncmp('siglay', currDimsNames, length(dimnames{unlimdimID + 1}));
    nIdx = strncmp('node', currDimsNames, length(dimnames{unlimdimID + 1}));
    nt = currDimsLengths(tIdx);
    ns = currDimsLengths(sIdx);
    nd = currDimsLengths(nIdx);

    if strcmpi(varname, 'temp') || strcmpi(varname, 'salinity')
        % To make the scaling go from the initial value to the POLCOMS value,
        % we need to take the scale the difference between the end members by
        % the scaling factor at each time and add to the current time's value.
        sfvdata = nan(nd, ns, nt);
        ss = 0:1 / (nt - 1):1; % scale from 0 to 1.
        startdata = squeeze(data(:, :, 1)); % use the first modelled time step
        for tt = 1:nt;
            if tt == 1
                sfvdata(:, :, 1) = startdata;
            else
                td = fvtemp - startdata;
                sfvdata(:, :, tt) = startdata + (ss(tt) .* td);
            end
        end

        % Replace the values with the scaled interpolated values.
        netcdf.putVar(ncout, varid, sfvdata)
    else
        % We need to check if the dimension is unlimited, and use a start
        % and end with netcdf.putVar if it is. This is largely because
        % MATLAB can't handle unlimited dimensions in the same way as it
        % does finite dimensions.
        if wasUnlimited < 0
            % We can just dump the entire data without specifying over what
            % indices.
            netcdf.putVar(ncout, varid, data);
        else
            % Use the dimension length we extracted above to output the
            % data with the valid unlimited dimension format.
            netcdf.putVar(ncout, varid, zeros(length(currDimsLengths), 1), currDimsLengths, data);
        end
    end

    if ftbverbose
        fprintf('done.\n')
    end
end

netcdf.close(nc)
netcdf.close(ncout)

if ftbverbose
    fprintf('%s : export complete.\n', subname)
end

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end