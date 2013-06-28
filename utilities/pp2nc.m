function pp2nc(file, convsh, pp2nc_tcl)
% Child function to call the convsh program to convert the obscure pp
% format to a sensible NetCDF which we can more easily read.
% 
% This requires the convsh and xconv functions available from:
% 
%   http://badc.nerc.ac.uk/help/software/xconv/index.html
% 
% Follow the installation instructions for your platform (Linux, Windows,
% Mac etc.) prior to running this function.
% 
% This requires the pp2nc.tcl script in the utilities subdirectory of the
% MATLAB fvcom-toolbox.
% 
% INPUT:
%   file - cell array of PP file(s) to convert to NetCDF.
%   convsh - the path to the convsh binary.
%   pp2nc_tcl - the path to the TCL script (in the utilities fvcom-toolbox
%       file).
% 
% OUTPUT:
%   NetCDF files in the same directory as the input PP files but with a .nc
%   file extension.
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Karen Amoudry (National Oceanography Centre, Liverpool)
% 
% PWC Revision history:
%   2013-06-24 Extracted version from the get_MetUM_forcing.m script and
%   set as standalone version.
%
% KJA Revision history:
%   2013-06-25 Added support for paths with spaces in. Also added support
%   for creation of NetCDF files in Windows (convsh will only take a .nc
%   filename, not a whole path).

nf = length(file);

for ff = 1:nf
    if ~isnan(file{ff})
        if exist(file{ff}, 'file') ~= 2
            error('File %s not found', file)
        end

        [path, name, ~] = fileparts(file{ff});
        out = [name, '.nc'];

        goback = pwd;
        cd(path)
        % Warn if it failed for some reason
        [res, msg] = system([convsh, ' "', pp2nc_tcl, '" -i "', file{ff}, '" -o "', out, '"']);
        cd(goback)
        if res ~= 0
            warning('Converion of %s to NetCDF failed. Conversion output:\n%s', file{ff}, msg)
        end
    end
end
