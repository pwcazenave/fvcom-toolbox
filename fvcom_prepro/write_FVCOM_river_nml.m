function write_FVCOM_river_nml(Mobj, nml_file, nc_file, vString)
% Write a namelist for the river nodes.
%
% write_FVCOM_river_nml(Mobj, nml_file, nc_file)
%
% DESCRIPTION:
%   Using the output of get_POLCOMS_rivers in Mobj, output the river name,
%   input file, grid node and vertical distribution into a namelist (.nml
%   file).
% 
%   NOTE: if the vertical distribution is uniform, the string 'uniform'
%   will be output. If a non-uniform distribution is detected, then the
%   name list will need to be edited to reflect the distribution specified
%   in the sigma.dat file.
%
% INPUT:
%   Mobj - MATLAB mesh object with the river data.
%   nml_file - full path to the output namelist file.
%   nc_file - full path to the NetCDF file containing the river data.
%   vString - optional, pass a string (e.g. 'uniform') to write as the
%   RIVER_VERTICAL_DISTRIBUTION in the namelist, bypassing the automated
%   string generation.
% 
% OUTPUT:
%   Namelist for inclusion in the main FVCOM namelist (RIVER_INFO_FILE).
% 
% EXAMPLE USAGE:
%   write_FVCOM_river_nml(Mobj, 'casename_river.nml', 'casename_river.nc')
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2013-03-21 - First version.
%   2013-08-16 - Add optional fourth argument of a string to supply as the
%   RIVER_VERTICAL_DISTRIBUTION (e.g. 'uniform').
%   2013-10-16 - Fix the handling of the optional vString argument.
%
%==========================================================================

subname = 'write_FVCOM_river_nml';

global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

nr = length(Mobj.river_nodes);

f = fopen(nml_file, 'w');
assert(f >= 0, 'Error writing to %s. Check permissions and try again.', nml_file)

% Build the vertical distribution string. Round to 15 decimal places so the
% unique check works (hopefully no one needs that many vertical layers...).
vDist = roundn(abs(diff(Mobj.siglev)), -15);
if length(unique(vDist)) == 1
    vString = '''uniform''';
elseif nargin <= 3
    vString = char();
    for ii = 1:length(vDist)
        vString = [vString, sprintf('%f ', vDist(ii))];
    end
end

for r = 1:nr
    fprintf(f, ' &NML_RIVER\n');
    fprintf(f, '  RIVER_NAME          = ''%s'',\n', Mobj.river_names{r});
    fprintf(f, '  RIVER_FILE          = ''%s'',\n', nc_file);
    fprintf(f, '  RIVER_GRID_LOCATION = %d,\n', Mobj.river_nodes(r));
    fprintf(f, '  RIVER_VERTICAL_DISTRIBUTION = %s\n', vString);
    fprintf(f, '  /\n');
end

fclose(f);

if ftbverbose
    fprintf('end   : %s\n', subname)
end
