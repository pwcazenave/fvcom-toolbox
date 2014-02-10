function write_FVCOM_probes(nml_file, interval, probes)
% Write a namelist for the timeseries extraction at specifc nodes.
%
% write_FVCOM_probes(Mobj, nml_file, interval, varlist, locations)
%
% DESCRIPTION:
%   Write a namelist of variables to extract at given positions with the
%   interval specified.
%
% INPUT:
%   nml_file - full path to the output namelist file.
%   interval - output interval (in seconds)
%   probe - struct of structs with fields whose names are the probe title.
%   Each probe struct must have the following fields:
%       node        - unstructured grid node number.
%       file        - file name for the current probe's output. If the
%           'variable' field is a cell array of multiple variables, each
%           output file name will have the variable name appended (e.g.
%           'Newlyn.dat' becomes 'Newlyn_el.dat' for elevation time
%           series).
%       description - description of this output. If multiple variable are
%           specified in variable, the length of the description must match
%           the number of variables.
%       variable    - variable name to extract from FVCOM. Define as a cell
%           array for multiple variables. Output file names in the namelist
%           will be appended with the variable name. Variable names are the
%           internal FVCOM variable names (which is not always the same as
%           the output file variable name). Check mod_probe.F in subroutine
%           PROBE_STORE for the full list of supported variables. I've
%           reproduced that list in the NOTES below.
%       longname    - long name for the current variable. Must also match
%       the number of variables.
%
% OUTPUT:
%   Namelist for inclusion in the main FVCOM namelist (PROBES_FILE).
%
% EXAMPLE USAGE:
%   probes.newlyn.node = 1045;
%   probes.newlyn.file = 'newlyn_elev.dat';
%   probes.newlyn.description = 'Elevation at Newlyn';
%   probes.newlyn.variable = 'el';
%   probes.newlyn.longname = 'Surface elevation (m)';
%   write_FVCOM_river_nml(Mobj, 'casename_probes.nml', 600, probes)
%
% NOTES:
%   Available output variables are currently (FVCOM v3.1.6) limited to:
%       - el - surface elevation
%       - t1, s1 - temperature and salinity
%       - rho1 - density
%       - u, v, ua, va - velocity and depth-averaged velocities
%       - ww - vertical water velocity
%       - w - vertical velocity in sigma system
%       - q2, l, q2l - turbulent kinetic energy (TKE), TKE lengthscale and
%           TKE times TKE lengthscale
%       - km, kh, kq - turbulent eddy viscosity for momentum, scalars and
%           q2/q2l
%       - aice, vice, uice2, vice2 - ice parameters
%       - csed - suspended sediment concentration
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2014-02-10 - First version.
%
%==========================================================================

subname = 'write_FVCOM_probes';

global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

fnames = fieldnames(probes);
np = length(fnames);

f = fopen(nml_file, 'w');
assert(f >= 0, 'Error writing to %s. Check permissions and try again.', nml_file)

for p = 1:np
    if iscell(probes.(fnames{p}).variable)
        nv = length(probes.(fnames{p}).variable);
        assert(length(probes.(fnames{p}).description) == nv, 'Specify one description per output variable.')
        assert(length(probes.(fnames{p}).description) == nv, 'Specify one longname per output variable.')
    else
        nv = 1;
        vname = probes.(fnames{p}).variable;
        lname = probes.(fnames{p}).longname;
        desc = probes.(fnames{p}).description;
    end
    for v = 1:nv
        if nv > 1
            vname = probes.(fnames{p}).variable{v};
            lname = probes.(fnames{p}).longname{v};
            desc = probes.(fnames{p}).description{v};
        end
        [pathstr, name, ext] = fileparts(probes.(fnames{p}).file);
        name = sprintf('%s_%s%s', name, vname, ext);
        file = fullfile(pathstr, name);
        fprintf(f, '&NML_PROBE\n');
        fprintf(f, ' PROBE_INTERVAL = ''seconds=%.1f'',\n', interval);
        fprintf(f, ' PROBE_LOCATION = %i,\n', probes.(fnames{p}).node);
        fprintf(f, ' PROBE_TITLE = ''%s'',\n', file);
        fprintf(f, ' PROBE_DESCRIPTION = ''%s'',\n', desc);
        fprintf(f, ' PROBE_VARIABLE = ''%s'',\n', vname);
        fprintf(f, ' PROBE_VAR_NAME = ''%s''\n', lname);
        fprintf(f, '/\n');
    end
end

fclose(f);

if ftbverbose
    fprintf('end   : %s\n', subname)
end
