function pp2nc_subset(file, convsh, pp2nc_tcl, xextents, yextents)
% Child function to call the convsh program to convert the PP format to a
% sensible NetCDF which we can more easily read.
%
% This requires the convsh and xconv functions available from:
%
%   http://badc.nerc.ac.uk/help/software/xconv/index.html
%
% Follow the installation instructions for your platform (Linux, Windows,
% Mac etc.) prior to running this function.
%
% This requires the subset.tcl script in the utilities subdirectory of the
% MATLAB fvcom-toolbox.
%
% INPUT:
%   file - cell array of PP file(s) to convert to NetCDF. Include
%       sub-directory in path, can include full path name or path relative
%       to working directory.
%   convsh - the path to the convsh binary (NB this can't contain spaces).
%   pp2nc_tcl - the path to the TCL script (in the utilities fvcom-toolbox
%       file).
%   xextents - [startx, endx, xinc] = start longitude coordinate, end
%       coordinate and grid increment (> 0.1degrees).
%   yextents - [starty, endy, yinc] = start longitude coordinate, end
%       coordinate and grid increment (> 0.1degrees).
%
% OUTPUT:
%   NetCDF files in the same directory as the input PP files but with a .nc
%   file extension.
%
% USAGE (PWC): Convert the PP files to NetCDF.
%   files = {'/path/to/file1.pp', '/path/to/file2.pp', '/path/to/file3.pp'};
%   convsh = '/users/modellers/pica/bin/convsh';
%   tcl = '/users/modellers/pica/MATLAB/fvcom-toolbox/utilities/pp2nc.tcl';
%   pp2nc(files, convsh, tcl);
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Karen Amoudry (National Oceanography Centre, Liverpool)
%    Judith Wolf (National Oceanography Centre, Liverpool)
%
% PWC Revision history:
%   2013-06-24 Extracted version from the get_MetUM_forcing.m script and
%   set as standalone version.
%
% KJA Revision history:
%   2013-06-25 Added support for paths with spaces in. Also added support
%   for creation of NetCDF files in Windows (convsh will only take a .nc
%   filename, not a whole path).
%
% JW Revision history:
%   2013-08-13 Problem for paths with spaces in, especially using Windows
%   Add backslash after pwd to define directory
%   Make sure convsh is pointing to right path (cannot use paths with
%   spaces in command, use relative paths from starting directory,
%   put output file in same working directory

nf = length(file);

for ff = 1:nf
    if ~isnan(file{ff})
        if exist(file{ff}, 'file') ~= 2
            error('File %s not found', file)
        end

        [pathname, filename, ~] = fileparts(file{ff});
        out = [filename, '.nc'];
%        infile =  [filename, '.pp'];    % JW added clear definition of input file
%        goback = pwd;
%        goback = strcat(pwd,'\')       % JW - add backslash (for Windows),
%        to define directory
       cd(pathname)
        % Warn if it failed for some reason
        [res, msg] = system([...
            convsh, ' "', pp2nc_tcl, '" -i "', file{ff}, '" -o "', out, ...
            '" -xs ', num2str(xextents(1)), ...
            ' -xe ', num2str(xextents(2)), ...
            ' -xi ', num2str(xextents(3)), ...
            ' -ys ', num2str(yextents(1)), ...
            ' -ye ', num2str(yextents(2)), ...
            ' -yi ',num2str( yextents(3))...
            ]);
%        [res, msg] = system([convsh, ' "', pp2nc_tcl, '" -i "', infile, '" -o "', out, '"']);
       cd(goback)
        if res ~= 0
            warning('Conversion of %s to NetCDF failed. Conversion output:\n%s', file{ff}, msg)
        end
    end
end
