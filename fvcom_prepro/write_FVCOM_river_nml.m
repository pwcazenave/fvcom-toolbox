function write_FVCOM_river_nml(Mobj, nml_file, nc_file)
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
%   2013-08-15 - Removed lines adding "uniform" river forcing as this
%   doens't seem to work. (ROM)
%
%==========================================================================

subname = 'write_FVCOM_river_nml';

global ftbverbose;
if ftbverbose
    fprintf(['\nbegin : ' subname '\n'])
end

nr = length(Mobj.river_nodes);

f = fopen(nml_file, 'w');
if f < 0
    error('Error writing to %s. Check permissions and try again.', nml_file)
end

% Build the vertical distribution string. Round to 15 decimal places so the
% unique check works (hopefully no one needs that many vertical layers...).
vDist = roundn(abs(diff(Mobj.siglev)), -15);
% if length(unique(vDist)) == 1
%     vString = '''uniform''';
% else
    vString = char();
    for ii = 1:length(vDist)
        vString = [vString, sprintf('%f ', vDist(ii))];
    end
% end

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
