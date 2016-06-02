% Unit test for write_FVCOM_tsobc.
%
% DESCRIPTION:
%   Currently checks against a reference data set for the following:
%       - number of nodes in the output
%       - number of sigma layers in the output
%       - number of time steps in the output
%       - range of values in the node arrays
%
% It uses a simplified POLCOMS NetCDF file from January, 2001 as the base
% input. The mesh object (Mobj) contains the required input for
% get_POLCOMS_tsobc as well as a set of 'known good' results
% (Mobj.temperature, Mobj.salt and Mobj.ts_times) for comparison against
% the new result.
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2013-06-02 First version.
%
%==========================================================================

matlabrc
close all
clc

% Set up our test environment.
[base, subname] = fileparts(mfilename('fullpath'));
cd(base)
%%
cd /users/modellers/pica/Code/fvcom-toolbox/tests/fvcom_prepro/
base = './';

addpath(fullfile(base, '../../fvcom_prepro'))
addpath(fullfile(base, '../../utilities'))


% Make some synthetic data to user for the test so we can easily compare
% the results.

% A boring grid.
nt = 30; % about a month of daily data
nx = 1024; % initial guess for the number of nodes
nz = 21; % vertical levels (= layers + 1)
[Mobj.lon, Mobj.lat] = meshgrid(1:ceil(sqrt(nx)), 1:ceil(sqrt(nx)));
% Offset each row by half a grid width.
Mobj.lon(1:2:end) = Mobj.lon(1:2:end) - max(max(diff(Mobj.lon, [], 2) / 2));
% Add the missing corners.
Mobj.lon = [Mobj.lon; min(Mobj.lon); max(Mobj.lon)];
Mobj.lat = [Mobj.lat; max(Mobj.lat); min(Mobj.lat)];
Mobj.lon = Mobj.lon(:) + 10;
Mobj.lat = Mobj.lat(:) + 50;
dt = delaunayTriangulation(Mobj.lon, Mobj.lat);
Mobj.tri = dt.ConnectivityList;
% Update the node count.
nx = length(Mobj.lon);
Mobj.h = randn(nx, 1) * 50;
% Set the "west" and "east" nodes as the open boundaries.
nodes = 1:nx;
Mobj.nObs = 2;
Mobj.read_obc_nodes{1} = nodes(Mobj.lon == min(Mobj.lon));
Mobj.read_obc_nodes{2} = nodes(Mobj.lon == max(Mobj.lon));

% Number of open boundary nodes.
nobn = length([Mobj.read_obc_nodes{:}]);

% Have a look see.
clf
trimesh(Mobj.tri, Mobj.lon, Mobj.lat, 'color', 'k')
hold on
plot(Mobj.lon([Mobj.read_obc_nodes{:}]), Mobj.lat([Mobj.read_obc_nodes{:}]), 'r.')

% Make the vertical grid (nice and simple at first).
Mobj.siglev = repmat(sigma_geo(nz, 1), [nx, 1]);
Mobj.siglay = nan(nx, nz - 1);
for i = 1:nz - 1
    Mobj.siglay(:, i) = mean(Mobj.siglev(:, i:i+1), 2);
end

% Some times.
Mobj.ts_times = linspace(51400, 51400 + nt - 1, nt);

% Data to write to file. Warm and fresh on top, cold and salty below.
in_temp = repmat(15, nobn, nz - 1, nt);
in_temp(:, 1, :) = 18;
in_salt = repmat(35, nobn, nz - 1, nt);
in_salt(:, 1, :) = 33;

write_FVCOM_tsobc(subname, Mobj.ts_times, nz - 1, in_temp, in_salt, Mobj)

%% Read in the written and compare against the inputs.

outvars = {'obc_h', 'obc_temp', 'obc_salinity', 'siglev', 'siglay', 'obc_nodes'};
obc_nodes = [Mobj.read_obc_nodes{:}];
inpvars = {Mobj.h(obc_nodes), in_temp, in_salt, ...
    Mobj.siglev(obc_nodes, :), ...
    Mobj.siglay(obc_nodes, :), ...
    obc_nodes};
assert(length(outvars) == length(inpvars), 'Inconsistent length of variables to compare.')
for v = 1:length(outvars)
    dump = double(ncread(sprintf('%s_tsobc.nc', subname), outvars{v}));
    % Try transposing vectors, everything else needs to be the same shape.
    if isvector(dump)
        try
            ncdiffs = dump - inpvars{v};
        catch
            ncdiffs = dump - inpvars{v}';
        end
    else
        ncdiffs = dump - inpvars{v};
    end
    if max(ncdiffs(:)) ~= 0
        warning('Input and output of ''%s'' differ by at most %g. Test FAILED.', outvars{v}, max(ncdiffs(:)))
    end
end
