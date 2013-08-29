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
        [varname, ~, ~, ~] = netcdf.inqVar(nc, ii - 1);
        
        if ismember(varname, varlist) || nargin == 1
            varid = netcdf.inqVarID(nc, varname);
            MetUM.(varname) = netcdf.getVar(nc, varid);
        end
    end
    
    if ftbverbose
        fprintf('done.\n')
    end

end

if ftbverbose
    fprintf('end   : %s \n', subname)
end