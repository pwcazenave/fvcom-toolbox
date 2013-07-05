function ncdata = get_POLCOMS_netCDF(files, varlist)
% Read temperature and salinity from NetCDF model output files and
% interpolate onto the open boundaries in Mobj.
%
% function struct = get_POLCOMS_netCDF(Mobj, files, varlist)
%
% DESCRIPTION:
%    Extract variables in varlist to a struct from the files given in
%    files.
%
% INPUT:
%   files   = Cell array of NetCDF file(s).
%   varlist = Cell array of variables names.
%
% OUTPUT:
%    Struct in which the field names are the variable names supplied in var
%    list. Each field has a data and units array which contain the data and
%    units variables from the NetCDF. If no units are found, it is left
%    blank for that variable.
%
% EXAMPLE USAGE
%    S = get_POLCOMS_netCDF({'/tmp/2000.nc', '/tmp/2001.nc', {'temp', 'salt'})
%
% NOTES:
%
%   - If you supply multiple files, there are a few assumptions:
%
%       - Variables are only appended if there are 3 or 4 dimensions; fewer
%       than that, and the values are assumed to be static across all the
%       given files (e.g. longitude, latitude etc.). The last dimension
%       is assumed to be time.
%       - The order of the files given should be chronological.
% 
%   - This has been tested on NetCDF files generated from the PML
%   POLCOMS-ERSEM daily mean model output
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2013-02-08 First version based on restart_FVCOM_AMM.m and
%    get_AMM_tsobc.m.
%
%==========================================================================

subname = 'get_POLCOMS_netCDF';

global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

% Get the results. Check we have a cell array, and if we don't, assume it's
% a file name.
if iscell(files)
    todo = length(files);
else
    todo = 1;
end

for ii = 1:todo

    if iscell(files)
        ftn = files{ii};
    else
        ftn = files;
    end

    if ftbverbose
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
        varid = netcdf.inqVarID(nc, getVar);

        data = double(netcdf.getVar(nc, varid, 'single'));
        if ii == 1
            ncdata.(getVar).data = data;
        else
            if ndims(data) < 3
                if strcmpi(getVar, 'time')
                    % If the dimension is time, we need to be a bit more
                    % clever since we'll need a concatenated time series
                    % (in which values are continuous and from which we
                    % can extract a sensible time). As such, we need to add
                    % the maximum of the existing time. On the first
                    % iteration, we should save ourselves the base time
                    % (from the units of the time variable).
                    ncdata.(getVar).data = [ncdata.(getVar).data; data + max(ncdata.(getVar).data)];
                else
                    % This should be a fixed set of values (probably lon or
                    % lat) in which case we don't need to append them, so
                    % just replace the existing values with those in the
                    % current NetCDF file.
                    ncdata.(getVar).data = data;
                end
            elseif ndims(data) == 3 || ndims(data) == 4
                % Concatenate along the last dimension and hope/assume it's
                % a time dimension.
                ncdata.(getVar).data = cat(ndims(data), ncdata.(getVar).data, data);
            else
                error('Unsupported number of dimensions in PML POLCOMS-ERSEM data')
            end
        end
        % Try to get some units (important for the calculation of MJD).
        try
            if ii == 1
                units = netcdf.getAtt(nc, varid, 'units');
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
        ncdata.(getVar).units = units;
    end

    netcdf.close(nc)

    if ftbverbose
        fprintf('done.\n')
    end

end

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end
