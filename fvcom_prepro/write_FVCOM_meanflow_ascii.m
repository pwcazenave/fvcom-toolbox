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

%% _meanflow.dat
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

%% _tide_node.dat -- nodes along the open boundaries.
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

%% _tide_cell.dat -- elements which have two nodes on an open boundary.
f = fopen([casename, '_tide_cell.dat'], 'w');
if f < 0
    error('Problem writing to .dat file. Check permissions and try again.')
end
if ~isfield(Mobj, 'read_obc_elements')
    error('Missing list of boundary element IDs. Run find_boundary_elements and try again.')
end
% Boundary element IDs
ne = Mobj.nObcElements;
fprintf(f, '%i\n', ne);
for j = 1:length(Mobj.read_obc_nodes); % number of boundaries
    for i = 1:numel(Mobj.read_obc_elements{j})
        fprintf(f, '%8i\n', Mobj.read_obc_elements{j}(i));
    end
end

fclose(f);

%% _tide_el.dat
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

%% _tide_uv.dat -- boundary velocities
% The format here is pretty funky. According to Dima's wrf_elj_obc.m
% script, for each time step, there's lines of depth averaged u and v
% followed by all the u components at each vertical level, then all the v
% components at each vertical level. All lines are prefixed with the
% current time in seconds relative to the start of the model (or mean
% flow? God knows).

f = fopen([casename, '_tide_uv.dat'], 'w');
if f < 0
    error('Problem writing to .dat file. Check permissions and try again.')
end
if ~isfield(Mobj, 'velocity')
    error('Missing mean flow velocity. Run get_POLCOMS_meanflow and try again.')
end

% Number of elements in the boundaries.
ne = Mobj.nObcElements;

% Number of time steps.
[~, nt] = size(Mobj.surfaceElevation);

% Do the depth averaged u then v for all nodes prefixed by the current
% time. So, wrap the whole shebang in a loop through time.
for t = 1:nt
    
    % Create a format string for the current time plus the number of
    % boundary elements.
    s = '%i';
    for ss = 1:ne
        if ss < ne
            s = [s, '%.4f\t'];
        else
            s = [s, '%.4f\n'];
        end
    end

    % Dump the time and mean u and then mean v vectors.
    fprintf(f, s, Mobj.meanflow_u(:, i));
    fprintf(f, s, Mobj.meanflow_v(:, i));

if ftbverbose
    fprintf('end   : %s \n', subname)
end


