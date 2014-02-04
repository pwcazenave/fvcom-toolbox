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
%       replaced. If the length of a given field is 1 (i.e. it's a single
%       value), that value will be repeated for the size of the fv_restart
%       variable of the same name. This is useful for writing constant
%       values to a spatially and/or temporally varying variable.
%
% OPTIONAL INPUTS:
%   The following options are argument-value pairs:
%       'new_times' - reset the restart file times to this time series
%           ([YYYY, MM, DD, HH, MM, SS] * number of timesteps).
%       'crop' - true/false switch. If true, then the restart
%       file data is cropped to the new_times time series; a value of false
%       ramps from the input example restart variable values. Defaults to
%       true.
%
% OUTPUT:
%   FVCOM restart file.
% 
% EXAMPLE USAGE:
%   Replace temperature and salinity:
%       indata.temp = interpolated_temp;
%       indata.salinity = interpolated_salinity;
%       write_FVCOM_restart('/tmp/fvcom_restart.nc', ...
%           '/tmp/fvcom_restart_interp.nc', indata)
%
%   Replace temperature and salinity and the times:
%       indata.temp = interpolated_temp;
%       indata.salinity = interpolated_salinity;
%       write_FVCOM_restart('/tmp/fvcom_restart.nc', ...
%           '/tmp/fvcom_restart_interp.nc', indata, ...
%           'new_times', inputConf.startDate)
%
%   Replace all salinity values with zeros:
%       indata.salinity = 0;
%       write_FVCOM_restart('/tmp/fvcom_restart.nc', ...
%           '/tmp/fvcom_restart_interp.nc', indata)
% 
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
% 
% Revision history:
%   2013-02-08 First version.
%   2013-02-15 Fix bug wherein only the last field in the new data would
%   only be added to the output netCDF file.
%   2013-03-13 Make the time rewriting optional and not just commented out.
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
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end

% Parse the input arguments. Set to new_times array to empty by default and
% only change if we're told otherwise in the input keyword-value pairs.
% Default to cropping the input in time rather than ramping up. This
% assumes you input data varying time. Mine (Pierre) don't (I only have a
% single snapshot which I want to be the final time step in the restart
% file).
crop = false;
new_times = [];
if nargin > 3
    for v = 1:2:length(varargin) - 1
        key = lower(varargin{v});
        val = varargin{v + 1};

        switch key
            case 'new_times'
                new_times = val;
            case 'crop'
                crop = val;
            otherwise
                warning('Unrecognised input argument keyword ''%s'' and value ''%s''.', key, varargin{v + 1})
        end
    end
end

% Set the dates output array empty unless we have different data in
% new_times (see loop below).
out_date = [];

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
    if ii ~= unlimdimID + 1 % netCDF indices start at zero
        dimid(ii) = netcdf.defDim(ncout, dimname, dimlen);
    else
        dimid(ii) = netcdf.defDim(ncout, dimname, netcdf.getConstant('NC_UNLIMITED'));
        
        % We're in time (since only time is unlimited), so sort out the
        % dates here. If we don't have new times, we'll check when we're
        % writing the time variables if out_dates is empty, and if it is,
        % just write the existing data.
        if ~isempty(new_times)
            % We need to process the time data a bit more carefully if we've
            % been asked to replace it. new_times can be a single date, in which
            % case we'll have to bracket it, or it can be an array of times, in
            % which case it must be the same length as the existing data.
            % Figure all that out here.
            if numel(new_times(:, 1)) == 1
                % Single time step, so bracket it by half an hour either way.
                out_date = [datenum(...
                        new_times(1), new_times(2), new_times(3), ...
                        new_times(4), new_times(5)-30, new_times(6)); ...
                    datenum(...
                        new_times(1), new_times(2), new_times(3), ...
                        new_times(4), new_times(5), new_times(6));...
                    datenum(...
                        new_times(1), new_times(2), new_times(3), ...
                        new_times(4), new_times(5)+30, new_times(6))
                    ];
            elseif numel(new_times(:, 1)) == dimlen
                out_date = new_times;
            else
                error('Replacement dates in ''new_times'' must be either a single date or an array of dates of the same size as the restart file duration.')
            end

            % Replace the dimlen with the number of new dates.
            dimlen = size(out_date, 1);
        end
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

    % Since the restart file has a number of time values, we'll ramp up the
    % replacement data from the existing start condition to the actual
    % value over the time steps. So, we need to know how many time steps we
    % actually have.

    % Get the dimension data ready for the replacement arrays.
    tIdx = strncmp(dimnames(unlimdimID + 1), currDimsNames, length(dimnames{unlimdimID + 1}));
    % Not sure about the hardcoded strings below...
    sIdx = strncmp('siglay', currDimsNames, length(dimnames{unlimdimID + 1}));
    nIdx = strncmp('node', currDimsNames, length(dimnames{unlimdimID + 1}));
    nt = currDimsLengths(tIdx);
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
            
            % Check the size of the new data. Assume that if the new data
            % is a scalar then we want to either:
            %   a) Replace the existing scalar with the new value (easy)
            %   b) Create an array of the same size as the existing restart
            %   file data, but with a value of the scalar.
            % a) is easy, b) is a little more compliated.
            if isscalar(data);
                % No need to ramp this if it's just a single value.
                sfvdata = indata.(fnames{vv});
            elseif isscalar(indata.(fnames{vv})) && ~isscalar(data)
                % Tile the scalar to the size of the input data.
                sfvdata = repmat(indata(fnames{vv}), size(data));

                % To make the scaling go from the initial value to the supplied
                % data value, we need to scale the difference between the end
                % members by the scaling factor at each time and add to the
                % current time's value.
            elseif crop
                % Limit the data to the time frame we've been given.
                sfvdata = indata.(fnames{vv})(:, :, time_idx);
            else
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
            end

            % Replace the values with the scaled interpolated values,
            % checking for unlimited dimensions as we go.
            if wasUnlimited < 0
                netcdf.putVar(ncout, varid, sfvdata)
            else
                netcdf.putVar(ncout, varid, zeros(length(currDimsLengths), 1), currDimsLengths, sfvdata)
            end

            writtenAlready = 1;

        % We might also want to replace the time. If so, supply the
        % keyword-value pairing of 'out_date' and a time array to replace
        % the times in the existing restart file with an arbitrary time
        % period.
        elseif strcmpi(varname, 'time') && writtenAlready == 0 && ~isempty(out_date)
            if ftbverbose
                fprintf('NEW DATA... ')
            end
            tmp_start_time = greg2mjulian(out_date(1), out_date(2), out_date(3) - 7, out_date(4), out_date(5), out_date(6));
            tmp_time = tmp_start_time:(tmp_start_time + nt - 1);
            netcdf.putVar(ncout, varid, tmp_time)

            writtenAlready = 1;

        elseif strcmpi(varname, 'Times') && writtenAlready == 0 && ~isempty(out_date)
            if ftbverbose
                fprintf('NEW DATA... ')
            end
            tmp_time = [];
            for i = 1:nt;
                tmp_time = [tmp_time, ...
                    sprintf('%-26s', ...
                    datestr(datenum(out_date(i, :)), ...
                        'yyyy-mm-dd HH:MM:SS.FFF'))];
            end
            % We have to specify the starting indices at the lengths of the
            % arrays because we've got an umlimited dimension.
            netcdf.putVar(ncout, varid, zeros(length(currDimsLengths), 1), currDimsLengths, tmp_time)

            writtenAlready = 1;

        elseif strcmpi(varname, 'Itime') && writtenAlready == 0 && ~isempty(out_date)
            if ftbverbose
                fprintf('NEW DATA... ')
            end
            tmp_start_time = greg2mjulian(out_date(1), out_date(2), out_date(3) - 7, out_date(4), out_date(5), out_date(6));
            tmp_time = tmp_start_time:(tmp_start_time + nt - 1);
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
        % We need to check if the dimension is unlimited, and use a
        % start and end with netcdf.putVar if it is. This is largely
        % because MATLAB can't handle unlimited dimensions in the same
        % way as it does finite dimensions.
        if wasUnlimited < 0
            % We can just dump the entire data without specifying over
            % what indices.
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
    fprintf(['end   : ' subname '\n'])
end
