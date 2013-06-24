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
% This uses the pp2nc.tcl script in the utilities subdirectory of the
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
% Revision history:
%   2013-06-24 Extracted version from the get_MetUM_forcing.m script and
%   set as standalone version.


if nargin == 1
    convsh = '/usr/local/bin/convsh';
end

nf = length(file);

for ff = 1:nf
    if ~isnan(file{ff})
        if exist(file{ff}, 'file') ~= 2
            error('File %s not found', file)
        end

        [path, name, ~] = fileparts(file{ff});
        out = fullfile(path, [name, '.nc']);

        system([convsh, ' ', pp2nc_tcl, ' -i ', file{ff}, ' -o ', out]);
    end
end
