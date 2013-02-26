function write_FVCOM_meanflow_ascii(Mobj, casename, data)
% Export mean flow forcing files hard-coded into mod_obcs2.F.
%
% function write_FVCOM_meanflow_ascii(Mobj, datfile, data)
%
% DESCRIPTION:
%    Setup an FVCOM hydrographic open boundary mean flow forcing file.
%
% INPUT:
%   Mobj     - MATLAB mesh object (with fields mf_time, siglay, siglev,
%               nObcNodes and read_obc_nodes).
%   casename - Output file prefix. Output files will be
%               /path/to/casename_suffix.dat, where suffix is one of:
%                   - meanflow
%                   - tide_cell
%                   - tide_node
%                   - tide_el
%                   - tide_uv
%   data     - 2D array of mean flow along the open boundary (nobc, time).
%
% OUTPUT:
%    FVCOM mean flow values along the FVCOM open boundary in a NETCDF file
%    named ncfile.
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2013-02-25 - First version.
% 
% TODO: Implement support for multiple open boundaries in all the outputs.
%
%==========================================================================

subname = 'write_FVCOM_meanflow_ascii';

global ftbverbose
if ftbverbose
    fprintf('\n'); fprintf(['begin : ' subname '\n']);
end

% _meanflow.dat
f = fopen([casename, '_meanflow.dat'], 'w');
if f < 0
    error('Problem writing to .dat file. Check permissions and try again.')
end
% Number of boundary nodes
fprintf(f, '%i\n', numel(Mobj.read_obc_nodes{1}));
% Boundary node IDs
for i = 1:numel(Mobj.read_obc_nodes{1})
    fprintf(f, '%i\n', Mobj.read_obc_nodes{1}(i));
end
% Sigma level distribution
s = '%i\t';
for ss = 1:length(Mobj.siglay)
    if ss < length(Mobj.siglay)
        s = [s, '%.4f\t'];
    else
        s = [s, '%.4f\n'];
    end
end
for i = 1:numel(Mobj.read_obc_nodes{1})
    fprintf(f, s, [i, abs(diff(Mobj.siglev))]);
end

% Number of times and boundary points
[nb, nt] = size(Mobj.velocity);

% Add the number of time steps
fprintf(f, '%i\n', nt);

s = '%.4f\n';
for ss = 1:nb
    if ss < nb
        s = [s, '%.4f\t'];
    else
        s = [s, '%.4f\n'];
    end
end
for i = 1:size(Mobj.velocity, 2)
    fprintf(f, s, [i - 1, Mobj.velocity(:, i)']);
end

fclose(f);

% _tide_cell.dat -- What's this? Element IDs?
f = fopen([casename, '_tide_cell.dat'], 'w');
if f < 0
    error('Problem writing to .dat file. Check permissions and try again.')
end
% Boundary node IDs
fprintf(f, '%i\n', numel(Mobj.read_obc_nodes{1}));
for i = 1:numel(Mobj.read_obc_nodes{1})
    fprintf(f, '%8i\n', Mobj.read_obc_nodes{1}(i));
end

fclose(f);

% _tide_node.dat
f = fopen([casename, '_tide_node.dat'], 'w');
if f < 0
    error('Problem writing to .dat file. Check permissions and try again.')
end
% Boundary node IDs
fprintf(f, '%i\n', numel(Mobj.obc_nodes(Mobj.obc_nodes ~= 0)));
for j = 1:length(Mobj.read_obc_nodes); % number of boundaries
    for i = 1:numel(Mobj.read_obc_nodes{j})
        fprintf(f, '%8i\n', Mobj.read_obc_nodes{j}(i));
    end
end

fclose(f);

% _tide_cell.dat
f = fopen([casename, '_tide_cell.dat'], 'w');
if f < 0
    error('Problem writing to .dat file. Check permissions and try again.')
end
if ~isfield(Mobj, 'read_obc_elements')
    error('Missing list of boundary element IDs. Run find_boundary_elements and try again.')
end
% Boundary element IDs
ne = 0;
for j = 1:length(Mobj.read_obc_elements)
    ne = ne + numel(Mobj.read_obc_elements{j});
end
fprintf(f, '%i\n', ne);
for j = 1:length(Mobj.read_obc_nodes); % number of boundaries
    for i = 1:numel(Mobj.read_obc_elements{j})
        fprintf(f, '%8i\n', Mobj.read_obc_elements{j}(i));
    end
end

fclose(f);

% _tide_el.dat
f = fopen([casename, '_tide_el.dat'], 'w');
if f < 0
    error('Problem writing to .dat file. Check permissions and try again.')
end
% Boundary node IDs
if ~isfield(Mobj, 'surfaceElevation')
    error('Missing predicted surface elevation necessary for mean flow.')
end
if ~isfield(Mobj, 'el_time')
    error('Missing predicted surface elevation time series necessary for mean flow.')
end

[nb, nt] = size(Mobj.surfaceElevation);

s = '%i\t';
for ss = 1:nb
    if ss < nb
        s = [s, '%.4f\t'];
    else
        s = [s, '%.4f\n'];
    end
end

for i = 1:nt
    fprintf(f, s', [round(Mobj.el_time(i)), Mobj.surfaceElevation(:, i)']);
end

fclose(f);

% _uv.dat -- boundary velocities? Are these tidal velocities rather than
% mean flow velocities, though? Hmmm...
f = fopen([casename, '_uv.dat'], 'w');
if f < 0
    error('Problem writing to .dat file. Check permissions and try again.')
end
if ~isfield(Mobj, 'velocity')
    error('Missing mean flow velocity. Run get_POLCOMS_meanflow and try again.')
end

% Number of elements in the boundaries.
ne = 0;
for j = 1:length(Mobj.read_obc_elements)
    ne = ne + numel(Mobj.read_obc_elements{j});
end
s = char();
for ss = 1:ne
    if ss < ne
        s = [s, '%.4f\t'];
    else
        s = [s, '%.4f\n'];
    end
end

for i = 1:ne
    fprintf(f, s, Mobj.meanflow_u(:, i));
    fprintf(f, s, Mobj.meanflow_v(:, i));
end

if ftbverbose
    fprintf('end   : %s \n', subname)
end


