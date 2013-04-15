function Mobj = get_POLCOMS_river_positions(Mobj, polcoms_grid, polcoms_ij)
% Parse the POLCOMS rivers data file.
%
% get_POLCOMS_rivers(Mobj, polcoms_file, polcoms_grid, polcoms_ij)
%
% DESCRIPTION:
%   Takes POLCOMS grid and index files and returns the positions and names
%   of the rivers within the index file.
%
% INPUT:
%   Mobj - MATLAB mesh object into which the outputs will be added.
%   polcoms_grid - NetCDF file of the POLCOMS grid.
%   polcoms_ij - indices in the POLCOMS grid at which each river is
%       located.
% 
% OUTPUT:
%   Mobj.river_xy - array of lon/lat positions of the POLCOMS rivers.
%   Mobj.river_polcoms_names - names for the rivers in the positions array.
%
% EXAMPLE USAGE:
%   Mobj = get_POLCOMS_rivers(Mobj, 'polcoms.nc', 'polcoms.index')
% 
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2013-04-15 - First version from the original get_POLCOMS_rivers.m.
%
%==========================================================================

subname = 'get_POLCOMS_river_positions';

global ftbverbose;
if ftbverbose
    fprintf(['\nbegin : ' subname '\n'])
end

% Check inputs
if exist(polcoms_grid, 'file') ~= 2
    error('file: %s does not exist or is not readable.', polcoms_grid)
end
if exist(polcoms_ij, 'file') ~= 2
    error('file: %s does not exist or is not readable.', polcoms_ij)
end

% Get the positions for each river from the NetCDF grid and the index file.
% Format is:
% 
% n
%   id1 i j Name1
%   id2 i j Name2
%   ...
%   idn i j Namen
fidx = fopen(polcoms_ij, 'r');
if fidx < 0
    error('Error opening index file (%s).', polcoms_ij);
end

c = 0; % line counter
while ~feof(fidx)
    line = fgetl(fidx);
    if isempty(line) || ~ischar(line)
        continue
    else
        c = c + 1;
    end
    
    if c == 1
        % First (valid) line should be number of rivers
        nridx = str2double(strtrim(line));
        % Preallocate the output arrays on the basis of this value. If the
        % file happens to have more lines, that shouldn't matter. Too few,
        % however, and you'll end up with NaNs at the end of your array. I
        % should probably remove them later...
        pc_idx = nan(nridx, 3);
        pc_name = cell(nridx, 1);
    else
        % We're in the data.
        S = regexpi(strtrim(line), ' +', 'split');
        pc_idx(c - 1, :) = [str2double(S{1}), str2double(S{2}), str2double(S{3})];
        pc_name{c - 1} = S{end};
    end
end
clear S line c

fclose(fidx);

% Now read in the NetCDF file and grab the real coordinates of those
% positions.
nc = netcdf.open(polcoms_grid, 'NOWRITE');
[~, numvars, ~, ~] = netcdf.inq(nc);

for ii = 1:numvars
    
    [varname, ~, ~, ~] = netcdf.inqVar(nc, ii - 1);
    varid = netcdf.inqVarID(nc, varname);

    if ftbverbose
        fprintf('\tvariable %s... ', varname)
    end
 
    switch varname
        case 'lon'
            pc_lon = netcdf.getVar(nc, varid);
            
        case 'lat'
            pc_lat = netcdf.getVar(nc, varid);

    end
    
    if ftbverbose
        fprintf('done.\n')
    end

end

clear numdims numvars dimnames 

netcdf.close(nc)

% Use the indices from the polcoms_ij file to get the real positions of the
% river nodes.
pc_riv_lonlat = [pc_lon(pc_idx(:, 2), 1), pc_lat(1, pc_idx(:, 3))'];

% Return the positions and names to the Mobj.
Mobj.rivers.positions = pc_riv_lonlat;
Mobj.rivers.names = pc_name;

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end
