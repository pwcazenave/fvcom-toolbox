function write_FVCOM_meanflow_ascii(Mobj, casename)
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
%            - Also requires two fields called meanflow_u and meanflow_v
%               which are arrays of u and v of sizes (nObcElements,
%               length(siglay), length(mf_times)).
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
    fprintf(['\nbegin : ' subname '\n'])
end

%% _meanflow.dat -- mean flow velocities at the open boundary elements (?).

% Create depth averaged velocity from the 3D velocity data in Mobj.
velocity = squeeze(mean(sqrt(Mobj.meanflow_u.^2 + Mobj.meanflow_v.^2), 2));

f = fopen([casename, '_meanflow.dat'], 'w');
if f < 0
    error('Problem writing to _meanflow.dat file. Check permissions and try again.')
end
% Number of boundary nodes
fprintf(f, '%8d\n', Mobj.nObcNodes);
% Boundary node IDs
for i = 1:Mobj.nObcNodes
    fprintf(f, '%8d\n', Mobj.read_obc_nodes{1}(i));
end
% Sigma level distribution
s = '%8d';
for ss = 1:length(Mobj.siglay)
    if ss < length(Mobj.siglay)
        s = [s, '%8.4f'];
    else
        s = [s, '%8.4f\n'];
    end
end
for i = 1:numel(Mobj.read_obc_nodes{1})
    fprintf(f, s, [i, abs(diff(Mobj.siglev))]);
end

% Number of times and boundary points
[nb, nt] = size(velocity);

% Add the number of time steps
fprintf(f, '%i\n', nt);

s = '%8.4f\n';
for ss = 1:nb
    if ss < nb
        s = [s, '%8.4f'];
    else
        s = [s, '%8.4f\n'];
    end
end
for i = 1:size(velocity, 2)
    fprintf(f, s, [i - 1, velocity(:, i)']);
end

fclose(f);

%% _tide_node.dat -- nodes along the open boundaries.
f = fopen([casename, '_tide_node.dat'], 'w');
if f < 0
    error('Problem writing to _tide_node.dat file. Check permissions and try again.')
end
% Boundary node IDs

% Get a list of the open boundary nodes. Transpose Mobj.obc_nodes so the
% order of the boundary nodes is preserved.
tmpObcNodes = Mobj.obc_nodes';
% Flip it back so it's the same shape as it would have been using the old
% code.
ObcNodes = tmpObcNodes(tmpObcNodes ~= 0)';

fprintf(f, '%8d\n', numel(ObcNodes));
for i = 1:numel(ObcNodes(i))
    fprintf(f, '%8i\n', ObcNodes(i));
end

fclose(f);

%% _tide_cell.dat -- elements which have two nodes on an open boundary.
f = fopen([casename, '_tide_cell.dat'], 'w');
if f < 0
    error('Problem writing to _tide_cell.dat file. Check permissions and try again.')
end
if ~isfield(Mobj, 'read_obc_elements')
    error('Missing list of boundary element IDs. Run find_boundary_elements and try again.')
end
% Boundary element IDs
ne = Mobj.nObcElements;
fprintf(f, '%8d\n', ne);
for j = 1:Mobj.nObs; % number of open boundaries
    for i = 1:numel(Mobj.read_obc_elements{j})
        fprintf(f, '%8i\n', Mobj.read_obc_elements{j}(i));
    end
end

fclose(f);

%% _tide_el.dat -- surface elevations with time.
f = fopen([casename, '_tide_el.dat'], 'w');
if f < 0
    error('Problem writing to _tide_el.dat file. Check permissions and try again.')
end
% Boundary node IDs
if ~isfield(Mobj, 'surfaceElevation')
    error('Missing predicted surface elevation necessary for mean flow.')
end
if ~isfield(Mobj, 'el_time')
    error('Missing predicted surface elevation time series necessary for mean flow.')
end

[nb, nt] = size(Mobj.surfaceElevation);

s = '%8d';
for ss = 1:nb
    if ss < nb
        s = [s, '%8.4f'];
    else
        s = [s, '%8.4f\n'];
    end
end

for i = 1:nt
    fprintf(f, s', [round(Mobj.el_time(i)), Mobj.surfaceElevation(:, i)']);
end

fclose(f);

%% _tide_uv.dat -- boundary velocities
% The format here is pretty funky. According to Dima's wrt_elj_obc.m
% script, for each time step, there's lines of depth averaged u and v
% followed by all the u and v components at each vertical level. All lines
% are prefixed with the current time in seconds relative to the start of
% the model (or mean flow time series? God knows).

f = fopen([casename, '_tide_uv.dat'], 'w');
if f < 0
    error('Problem writing to _tide_uv.dat file. Check permissions and try again.')
end

% Number of elements in the boundaries.
ne = Mobj.nObcElements;

% Number of time steps.
nt = length(Mobj.mf_times);

% Number of vertical layers.
nz = length(Mobj.siglay);

% Create a format string for the each time step plus the number of boundary
% elements.
s = '%8d';
for ss = 1:ne
    
    if ss < ne
        s = [s, '%8.4f'];
    else
        s = [s, '%8.4f\n'];
    end
end

% Do the depth averaged u then v for all nodes prefixed by the current
% time. So, wrap the whole shebang in a loop through time.
for t = 1:nt
    
    % Time since the start of the time series (in seconds).
    iint = (Mobj.mf_times(t) - Mobj.mf_times(1)) * 24 * 3600;
    
    % Dump the time and mean u and then mean v vectors.
    fprintf(f, s, [iint; mean(Mobj.meanflow_u(:, :, t), 2)]);
    fprintf(f, s, [iint; mean(Mobj.meanflow_v(:, :, t), 2)]);
    
    % Now, for each vertical layer, dump the u and v vectors, prefixed with
    % time.
    for zz = 1:nz
        fprintf(f, s, [iint; Mobj.meanflow_u(:, zz, t)]);
        fprintf(f, s, [iint; Mobj.meanflow_v(:, zz, t)]);
    end
end

fclose(f);

%% _elj_obc.dat -- surface elevation time series at open boundary nodes.

% This is almost identical to _tide_el.dat but lacks the time stamp in the
% first column.

f = fopen([casename, '_elj_obc.dat'], 'w');
if f < 0
    error('Problem writing to _elj_obc.dat file. Check permissions and try again.')
end

nt = size(Mobj.surfaceElevation, 2);

for t = 1:nt
    fprintf(f, '%8.4f', Mobj.surfaceElevation(:, t));
    fprintf(f, '\n');
end

if ftbverbose
    fprintf('end   : %s\n', subname)
end


