function MetUM = read_MetUM_forcing(files, varlist)
% Read Met Office Unified Model (TM) (hereafter MetUM) netCDF files and
% extact the variables in varlist.
%
% MetUM = read_MetUM_forcing(files, varlist)
%
% DESCRIPTION:
%   Given a cell array of file names, extract the variables in varlist into
%   a struct. Field names are the variable names gives.
%
% INPUT:
%   files - cell array of file names (full paths to the files)
%   varlist [optional] - list of variable names to extract. If omitted, all
%   variables are extracted.
%
% OUTPUT:
%   MetUM - MATLAB struct with the data from the variables specified in
%   varlist.
%
% EXAMPLE USAGE:
%   varlist = {'x', 'y', 'sh', 'x-wind', 'y-wind', 'field202'};
%   files = {'/tmp/metum/sn_2011010100_s00.nc', ...
%       '/tmp/metum/sn_201101016_s00.nc'};
%   MetUM = read_MetUM_forcing(files, varlist);
%
% TODO:
%   - Fix the generation of MetUM.t (currently some values are missing for
%   some reason.
%   - Find out why some of the variables have an inconsistent number of 3rd
%   dimension (vertical) levels. This shouldn't happen beceause the number
%   of vertical levels in a given variable shouldn't change across multiple
%   output files.
%   - Extract only the relevant vertical levels (this may fix the issue
%   above as we'll only have a single vertical level at that point).
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% 2013-08-29 First version.
%
%==========================================================================

subname = 'read_MetUM_forcing';

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

assert(iscell(files), 'List of files provided must be a cell array.')

MetUM = struct();

for f = 1:length(files)
    
    if ftbverbose
        fprintf('File %s... ', files{f})
    end

    nc = netcdf.open(files{f}, 'NOWRITE');
    
    % Query the netCDF file to file the variable names. If the name matches
    % one in the list we've been given (or if we haven't been given any
    % particular variables), save it in the output struct.
    [~, numvars, ~, ~] = netcdf.inq(nc);

    for ii = 1:numvars
        % Find the name of the current variable
        [varname, ~, ~, varAtts] = netcdf.inqVar(nc, ii - 1);

        if ismember(varname, varlist) || nargin == 1
            varid = netcdf.inqVarID(nc, varname);

            % Some variables contain illegal (in MATLAB) characters. Remove
            % them here.
            safename = regexprep(varname, '-', '');

            % Append the data on the assumption the last dimension is time.
            % Don't append data with only 2 dimensions as it's probably
            % longitude or latitude data. The time variable ('t') is
            % turned into a list of time stamps.
            tmpdata = squeeze(netcdf.getVar(nc, varid, 'double'));
            nn = ndims(tmpdata);

            if isfield(MetUM, safename)
                switch varname
                    case {'x', 'y', 'x_1', 'y_1'}
                        continue
                    case {'t', 't_1'}
                        % Get the time attribute so we can store proper
                        % times.
                        tt = fix_time(nc, varid, varAtts);
                        MetUM.(safename) = cat(1, MetUM.(safename), tt);
                    otherwise
                        try
                            % Append along last dimension
                            MetUM.(safename) = cat(nn, MetUM.(safename), tmpdata);
                        catch
                            fprintf('\n')
                            warning('Couldn''t append %s to the existing field from file %s.', safename, files{f})
                        end

                end
            else % first time around
                MetUM.(safename) = tmpdata;
                if strcmpi(varname, 't') || strcmpi(varname, 't_1')
                    % Get the time attribute so we can store proper times.
                    MetUM.(safename) = fix_time(nc, varid, varAtts);
                end
            end
        end
    end

    if ftbverbose
        fprintf('done.\n')
    end

end

% Squeeze out singleton dimensions.
fields = fieldnames(MetUM);
for i = 1:length(MetUM)
    MetUM.(fields{i}) = squeeze(MetUM.(fields{i}));
end

if ftbverbose
    fprintf('end   : %s \n', subname)
end

function tt = fix_time(nc, varid, varAtts)
% Little helper function to get the time attribute so we can store proper
% times.
%
% INPUT:
%   nc : netCDF file handle
%   varid : current variable ID
%   varAtts : number of variable attributes
%
% OUTPUT:
%   tt : time string for the current file
%
for j = 1:varAtts
    timeatt = netcdf.inqAttName(nc, varid, j - 1);
    if strcmpi(timeatt, 'time_origin')
        timeval = netcdf.getAtt(nc, varid, timeatt);
    end
end
mt = datenum(timeval, 'dd-mmm-yyyy:HH:MM:SS');
% t is the offset in days relative to mt.
t = netcdf.getVar(nc, varid, 'double');
if isempty(t)
    % This isn't right. Might be easier to get the time from the file
    % name...
    t = 0;
end
% Add the offset and convert back to a date string.
tt = datestr(mt + t, 'yyyy-mm-dd HH:MM:SS');
