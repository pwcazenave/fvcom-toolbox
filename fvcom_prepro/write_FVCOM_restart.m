function write_FVCOM_restart(fv_restart, out_restart, indata, varargin)
% Duplicate an FVCOM restart file, replacing variable values with those
% specified in the struct data. 
% 
% function write_FVCOM_restart(fv_restart, out_restart, indata)
% 
% DESCRIPTION:
%   Use an existing FVCOM restart file as a template, export all existing
%   data except for variables whose names match the data in the struct
%   'data'.
% 
% INPUT:
%   fv_restart  = full path to an existing FVCOM restart file.
%   out_restart = full path to the restart file to be created.
%   indata      = struct whose field names are the variable names to be
%   replaced.
% OPTIONAL INPUT (keyword-value pairs):
%   'out_date'  = [optional] reset the restart file times to this date
%   ([YYYY, MM, DD, HH, MM, SS]). This must be a single date only. If new
%   data are being provided, they must be a single time step only; existing
%   data will use the last time step in the restart file. The output file
%   will include three time steps to bracket the specified time by 30
%   minutes each way (to allow FVCOM some wiggle room when loading the
%   data). The data will be replicated over the three time steps.
%
% OUTPUT:
%   FVCOM restart file.
% 
% EXAMPLE USAGE:
%
%   Replace temperature and salinity but leave the times as is:
%       indata.temp = interpolated_temp;
%       indata.salinity = interpolated_salinity;
%       write_FVCOM_restart('/tmp/fvcom_restart.nc', ...
%           '/tmp/fvcom_restart_interp.nc', indata)
%
%   Replace temperature only and change the restart times to bracketed
%   values:
%       indata.temp = interpolated_temp;
%       write_FVCOM_restart('/tmp/fvcom_restart.nc', ...
%           '/tmp/fvcom_restart_interp.nc', indata, 'out_date', ...
%           [2003, 05, 25, 13, 34, 07])
%
%   Replace all temperatures with a single value leaving times as they are:
%       indata.temp = 13;
%       write_FVCOM_restart('/tmp/fvcom_restart.nc', ...
%           '/tmp/fvcom_restart_interp.nc', indata)
%
%   Replace temperature with a constant at all positions and times and all
%   times with a new array of times (for a restart file with five time
%   steps):
%       newtime = datevec(datenum([2004, 06, 20, 11, 09, 44]) + ...
%           [-0.5, -0.25, 0, 0.25, 0.5]);
%       indata.temp = 13;
%       write_FVCOM_restart('/tmp/fvcom_restart.nc', ...
%           '/tmp/fvcom_restart_interp.nc', indata, 'out_time', newtime)
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%   Karen Amoudtry (National Oceanography Centre, Liverpool)
%
% Revision history:
%   2013-02-08 First version.
%   2013-02-15 Fix bug wherein only the last field in the new data would
%   only be added to the output netCDF file.
%   2013-03-13 Make the time rewriting optional and not just commented out.
%   2014-02-04 Incorporate Karen's functionality (see revision history
%   below), but with the ability to retain the existing behaviour (where a
%   new start time is still optional). User-specified constants are also
%   supported but instead of specifying a new input argument, if a single
%   scalar value is given in the input struct but the output is non-scalar
%   (i.e. an array), then that scalar is tiled to the size of the expected
%   output array.
%
% KJA Revision history:
%   2014-01-23 Add functionality to specify length of time series in output
%   file.
%   2014-01-31 Add functionality to replace user-specified variables in the
%   output file with user-specified constant values.
%
%==========================================================================

subname = 'write_FVCOM_restart';

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

% Default to not bracketing data. If we've got a single time input, then
% this will be set to true and the data will be bracketed around that time.
% If we've got a new input time array, the new data will be ramped to the
% new input data values and bracketed is left false. If we have no new
% times, then the new data will also be ramped up.
bracketed = false;

if nargin > 3
    assert(rem(length(varargin), 2) == 0, 'Incorrect keyword-value arguments.')
    for aa = 1:2:length(varargin)
        key = varargin{aa};
        val = varargin{aa + 1};
        switch key
            case 'out_date'
                if isscalar(val(:, 1))
                    % Bracket the date by 30 minutes either way.
                    bracketed = true;
                    tOffset = 30;
                    out_date = val;
                    new_date = datevec([...
                        datenum(...
                            val(1), val(2), val(3), val(4), val(5) - tOffset, val(6)); ...
                        datenum(...
                            val(1), val(2), val(3), val(4), val(5), val(6));...
                        datenum(...
                        val(1), val(2), val(3), val(4), val(5) + tOffset, val(6))...
                        ]);
                else
                    % Assume we've been given a time series. We need to
                    % make sure (later) that this is the same length as the
                    % input file times.
                    new_date = val;
                end
            otherwise
                warning('Unrecognised input argument keyword ''%s'' and value ''%s''.', key, varargin{v + 1})
        end
    end
end

% Get the fieldnames which must match the variable names to be replaced
% (case sensitive).
fnames = fieldnames(indata);
nf = length(fnames);

nc = netcdf.open(fv_restart, 'NOWRITE');
ncout = netcdf.create(out_restart, 'clobber');

[numdims, numvars, numglobatts, unlimdimID] = netcdf.inq(nc);

% Define the dimensions for all variables.
dimid = nan(numdims, 1);
dimnames = cell(numdims, 1);
dimlengths = nan(numdims, 1);
for ii = 1:numdims
    [dimname, dimlen] = netcdf.inqDim(nc, ii - 1);
    % If we've been asked to rewrite the times, then set the length
    % of the time dimension to three (bracketed by two time steps each
    % way).
    if exist('out_date', 'var') && ii == unlimdimID + 1
        if length(out_date(:, 1)) == dimlen
            % We're replacing times with those in new_date
            out_date = new_date;
        elseif isvector(out_date)
            % We have bracketed dates
            out_date = new_date;
        end
        dimlen = length(out_date(:, 1));
    end

    if ii ~= unlimdimID + 1 % netCDF indices start at zero
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
    
    % Put that back into the output netCDF file.
    netcdf.putAtt(ncout, netcdf.getConstant('NC_GLOBAL'), gattname, gattval);
end

netcdf.endDef(ncout);

% Get the existing data and output to the new netCDF file, except for
% variables which match the fieldnames in the data struct.
for ii = 1:numvars

    [varname, ~, varDimIDs, ~] = netcdf.inqVar(nc, ii - 1);
    varid = netcdf.inqVarID(nc, varname);

    if ftbverbose
        fprintf('\tvariable %s... ', varname)
    end

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

    % Since the restart file has a number of time values, we'll ramp up the
    % replacement data from the existing start condition to the actual
    % value over the time steps. So, we need to know how many time steps we
    % actually have.

    % Get the dimension data ready for the replacement arrays.
    tIdx = strncmp(dimnames(unlimdimID + 1), currDimsNames, length(dimnames{unlimdimID + 1}));
    nt = currDimsLengths(tIdx);
    % Not sure about the hardcoded strings below...
    sIdx = strncmp('siglay', currDimsNames, length(dimnames{unlimdimID + 1}));
    nIdx = strncmp('node', currDimsNames, length(dimnames{unlimdimID + 1}));
    ns = currDimsLengths(sIdx);
    nd = currDimsLengths(nIdx);
    if isempty(nd)
       % We've got data on the elements (i.e. u and v)
       nIdx = strncmp('nele', currDimsNames, length(dimnames{unlimdimID + 1}));
       nd = currDimsLengths(nIdx);
    end
    
    % Iterate through the field names to see if we're on one of the
    % variables to be replaced.

    % Set variable so we know if we've already written this variable to the
    % output file.
    writtenAlready = 0;
    for vv = 1:nf
        if strcmp(varname, fnames{vv}) && writtenAlready == 0
            if ftbverbose
                fprintf('NEW DATA... ')
            end

            % Only repeat the values when we've bracketed the data,
            % otherwise scale up from the first time step.
            if bracketed
                % Use the new input data repeated the required number of
                % times.
                if wasUnlimited < 0 % no time, easy peasy.
                    sfvdata = indata.(fnames{vv});
                else
                    % Damn. We've got to filter based on the shape of the
                    % input data first (do we have a scalar input we need
                    % to repeat to the shape of the output?) and then on
                    % the basis of the expected output (is it 1D, 2D or
                    % 3D?).
                    if isscalar(indata.(fnames{vv}))
                        if isempty(ns) && isempty(nd)
                            sfvdata = indata.(fnames{vv});
                        elseif isempty(ns) || isempty(nd)
                            if ftbverbose
                                fprintf('tiling input scalar to non-scalar array... ')
                            end
                            sfvdata = repmat(indata.(fnames{vv}), [nd, nt]);
                        elseif ~isempty(ns) && ~isempty(nd)
                            if ftbverbose
                                fprintf('tiling input scalar to non-scalar array... ')
                            end
                            sfvdata = repmat(indata.(fnames{vv}), [nd, ns, nt]);
                        else
                            error('Hmmm... FVCOM doesn''t have 4D arrays...')
                        end
                    else
                        % We don't have scalar input so repeat the last
                        % time in the given input nt times.
                        if isempty(ns) && isempty(nd)
                            sfvdata = indata.(fnames{vv});
                        elseif isempty(ns) || isempty(nd)
                            sfvdata = repmat(indata.(fnames{vv})(:, :, end), [1, nt]);
                        elseif ~isempty(ns) && ~isempty(nd)
                            sfvdata = repmat(indata.(fnames{vv})(:, :, end), [1, 1, nt]);
                        else
                            error('Hmmm... FVCOM doesn''t have 4D arrays...')
                        end
                    end
                    if ftbverbose
                        fprintf('clipping in time... ')
                    end
                end
            else
                % To make the scaling go from the initial value to the
                % supplied data value, we need to scale the difference
                % between the end members by the scaling factor at each
                % time and add to the current time's value.
                data = netcdf.getVar(nc, varid);
                if ~isscalar(data)
                    if isscalar(indata.(fnames{vv})) && ftbverbose
                        fprintf('tiling input scalar to non-scalar array... ')
                    end
                    sfvdata = nan(nd, ns, nt);
                    ss = 0:1 / (nt - 1):1; % scale from 0 to 1.
                    startdata = squeeze(data(:, :, 1)); % use the first modelled time step
                    for tt = 1:nt
                        if tt == 1
                            sfvdata(:, :, 1) = startdata;
                        else
                            td = indata.(fnames{vv}) - startdata;
                            sfvdata(:, :, tt) = startdata + (ss(tt) .* td);
                        end
                    end
                    if ftbverbose
                        fprintf('ramping data in time... ')
                    end
                else
                    sfvdata = data;
                end
            end
            % Replace the values with the scaled interpolated values,
            % checking for unlimited dimensions as we go.
            if wasUnlimited < 0
                netcdf.putVar(ncout, varid, sfvdata)
            else
                netcdf.putVar(ncout, varid, zeros(length(currDimsLengths), 1), currDimsLengths, sfvdata)
            end

            writtenAlready = 1;

        % We might also want to replace the time. If so, supply a fourth
        % argument (out_date) to replace the times in the existing
        % restart file with an arbitrary time period.
        elseif strcmpi(varname, 'time') && writtenAlready == 0 && exist('out_date', 'var')
            if ftbverbose
                fprintf('NEW DATA... ')
            end
            tmp_time = greg2mjulian(out_date(:, 1), out_date(:, 2), out_date(:, 3), out_date(:, 4), out_date(:, 5), out_date(:, 6));
            netcdf.putVar(ncout, varid, tmp_time)

            writtenAlready = 1;

        elseif strcmpi(varname, 'Times') && writtenAlready == 0 && exist('out_date', 'var')
            if ftbverbose
                fprintf('NEW DATA... ')
            end
            tmp_time = [];
            for i = 1:nt;
                tmp_time = [tmp_time, sprintf('%-026s', datestr(datenum(out_date(i, :)), 'yyyy-mm-dd HH:MM:SS.FFF'))];
            end
            netcdf.putVar(ncout, varid, tmp_time)

            writtenAlready = 1;

        elseif strcmpi(varname, 'Itime') && writtenAlready == 0 && exist('out_date', 'var')
            if ftbverbose
                fprintf('NEW DATA... ')
            end
            tmp_time = greg2mjulian(out_date(:, 1), out_date(:, 2), out_date(:, 3), out_date(:, 4), out_date(:, 5), out_date(:, 6));
            netcdf.putVar(ncout, varid, floor(tmp_time))

            writtenAlready = 1;

        end
    end

    % If writtenAlready is zero, we haven't had one of the variables we're
    % replacing, so just dump the existing data.
    if writtenAlready == 0
        if ftbverbose
            fprintf('existing data... ')
        end

        % Grab the data.
        data = netcdf.getVar(nc, varid);

        % We need to check if the dimension is unlimited, and use a
        % start and end with netcdf.putVar if it is. This is largely
        % because MATLAB can't handle unlimited dimensions in the same
        % way as it does finite dimensions.
        if wasUnlimited < 0
            % We can just dump the entire data without specifying over
            % what indices.
            netcdf.putVar(ncout, varid, data);
        else
            % If we're clipping in time, use only the last value repeated
            % three times for the bracketed times. Otherwise, we just dump
            % it as is. FVCOM currently has no 4D variables: since the two
            % spatial dimensions are collapsed into one, the maximum number
            % of dimensions is 3: time, depth and position.
            if any(ismember(currDimsNames, 'time')) && exist('out_date', 'var')
                if ftbverbose
                    fprintf('clipping in time... ')
                end
                if isvector(data) % 1D data
                    data = data(end - (nt - 1):end);
                elseif ismatrix(data) % 2D data
                    data = repmat(data(:, end), [1, nt]);
                else % 3D data
                    data = repmat(data(:, :, end), [1, 1, nt]);
                end
            end

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
    fprintf('end   : %s \n', subname)
end
