% Create the necessary files for FVCOM from a given SMS input unstructured
% grid.
%
% This version uses the following forcings:
%   - E-HYPE river climatology discharge based on CEH rivers
%   - EA river temperature climatology
%   - River salinity based on my river morphology classification
%   - CFS surface forcing (wind, heat, preciptation/evaporation)
%   - HYCOM open boundary temperature/salinity and initial conditions
%   - TPXO predicted tidal elevations at the boundary
%
% This requires:
%   - my fork of the fvcom-toolbox
%       > https://gitlab.em.pml.ac.uk/pica/fvcom-toolbox
%   - the Tide Model Driver toolbox and data
%       > http://polaris.esr.org/ptm_index.html
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%
%   2015-05-15 Poached the Irish Sea domain version of this script. Removed
%   the version history as it was a bit long and pointless.
%   2015-08-14 Use ECMWF-ERA20C forcing instead of NCEP.

% matlabrc
close all
clc

cd('/users/modellers/rito/Models/git/fvcom/locate-implementation/loc_v0/grids')

global ftbverbose
ftbverbose = 1; % be noisy

addpath('/users/modellers/rito/Models/git/fvcom/fvcom-toolbox/utilities')
addpath('/users/modellers/rito/Models/git/fvcom/fvcom-toolbox/fvcom_prepro/')
% addpath('/users/modellers/pica/Code/MATLAB/func/')
addpath('/users/modellers/rito/matlab/air-sea')
conf.obc_tides.dirTMD = '/users/modellers/pica/Code/MATLAB/toolboxes/TMD/';
addpath(conf.obc_tides.dirTMD)
addpath(fullfile(conf.obc_tides.dirTMD, 'FUNCTIONS'))
addpath('scriptlets')

%%%------------------------------------------------------------------------
%%%                          INPUT CONFIGURATION
%%%
%%% Define the model input parameters here e.g. number of tidal components,
%%% type of boundary forcing, estimated current velocity, sponge layer
%%% coefficient and radius, forcing source etc.
%%%
%%%------------------------------------------------------------------------

conf.base = '/users/modellers/rito/Models/git/fvcom/locate-implementation/loc_v0/grids';

% Which version of FVCOM are we using (for the forcing file formats)?
conf.FVCOM_version = '3.2.2';


%%%------------------------------------------------------------------------
%%%                           Spatial stuff
%%%------------------------------------------------------------------------

% Case name for the model inputs and outputs
conf.casename = 'lyme_bay_v04';

conf.coordType = 'cartesian'; % 'cartesian' or 'spherical'
% Input grid UTM Zone (if applicable)
conf.utmZone = {'30 U'}; % syntax for utm2deg

% Give some names to the boundaries. This must match the number of node
% strings defined in SMS. Ideally, the order of the names should match the
% order in which the boundaries were made in SMS.
conf.boundaryNames = {'South'}; % Number is (relatively) important!

%%%------------------------------------------------------------------------
%%%                     Model constants and stuff
%%%------------------------------------------------------------------------

% Type of open boundary treatment (see Table 6.1 (?) in the manual).
% 1 - Active (ASL): sea level is specified at the open boundary.
% 2 - Clamped: zeta = 0 at the boundary (Beardsley and Haidvogel, 1981).
% 3 - Implicit Gravity Wave Radiation.
% 4 - Partial Clamped Gravity Wave Radiation (Blumberg and Kantha, 1985).
% 5 - Explicit Orlanski Radiation (Orlanski, 1976; Chapman, 1985)
conf.obc_type = 1;

% Sponge layer parameters
conf.sponge.radius = 15000; % in metres, or -1 for variable width.
conf.sponge.coeff = 0.001;

% z0 value in metres (for uniform) or 'variable' (requires a separate
% roughness (z0) distribution grid).
conf.bedRoughness = 0.035; % or 0.015, 0.025 or 0.03 - Davies and Jones (1990) shelf model

% Estimated velocity (m/s) and tidal range (m) for time step estimate
conf.estVel = 4;
conf.estRange = 14;

%%%------------------------------------------------------------------------
%%%                      END OF INPUT CONFIGURATION
%%%------------------------------------------------------------------------

%% Process the files needed for all months.

% Read the input mesh and bathymetry. Also creates the data necessary for
% the Coriolis correction in FVCOM.
Mobj = read_sms_mesh(...
    '2dm', fullfile(conf.base,conf.casename, [conf.casename, '.2dm']),...
    'coordinate', conf.coordType, 'addCoriolis', false);
%     'bath', fullfile(conf.base,conf.casename, [conf.casename, '.dat']),...

% Clean out some unneeded fields.
Mobj = rmfield(Mobj, 'riv_nodes');

% Convert the coordinates to whatever is missing.
if Mobj.have_lonlat == 0
    utmZones = cellfun(@(x) repmat(x, length(Mobj.x), 1), conf.utmZone, 'uni', false);
    [Mobj.lat, Mobj.lon] = utm2deg(Mobj.x, Mobj.y, utmZones{1});
    Mobj.have_lonlat = true;
    clear utmZones
end
if Mobj.have_xy == 0
    tmpZone = regexpi(conf.utmZone, '\ ', 'split');
    [Mobj.x, Mobj.y] = wgs2utm(Mobj.lon, Mobj.lon, str2double(char(tmpZone{1}(1))), 'N');
    Mobj.have_xy = true;
    clearvars tmpZone
end

% Add grid metrics.
Mobj = setup_metrics(Mobj);

% write ADMESH 14 format file


% Create a Coriolis file from the bathy which varies with latitude. Given
% the size of the domain, this is probably necessary.
if ~Mobj.have_cor
    Mobj = add_coriolis(Mobj, 'uselatitude');
end

% Parse the open boundary nodes and add accordingly.
if Mobj.have_strings
    for i = 1:size(Mobj.read_obc_nodes, 2)
        nodeList = double(cell2mat(Mobj.read_obc_nodes(i)));
        Mobj = add_obc_nodes_list(Mobj, nodeList, conf.boundaryNames{i}, conf.obc_type);
    end
    clear nodeList
end

% Create a sponge layer
if Mobj.have_strings
    for i = 1:size(Mobj.read_obc_nodes, 2)
        % Get the list of nodes in this boundary
        nodeList = double(cell2mat(Mobj.read_obc_nodes(i)));

        if conf.sponge.radius < 0 % if we want a variable sponge radius
            if i==1
                % Create an array to store the radii
                Mobj.sponge_rad = zeros(size(Mobj.sponge_nodes));
            end
            % calculate the sponge radius
            conf.spongeRadius = calc_sponge_radius(Mobj, nodeList);
            % Add the sponge nodes to the list
            Mobj = add_sponge_nodes_list(Mobj,nodeList,...
                [conf.boundaryNames{i}, ' sponge'], conf.sponge.radius,...
                conf.sponge.coeff);
        else
            Mobj = add_sponge_nodes_list(Mobj,nodeList,...
                [conf.boundaryNames{i}, ' sponge'], conf.sponge.radius,...
                conf.sponge.coeff);
        end
        clear nodeList
    end
end
[Mobj] = write_admesh_mesh(Mobj,'output_directory',fullfile(conf.base,conf.casename),'native_coord','spherical');
%  [Mobj] = write_admesh_mesh(Mobj,varargin)
save locate_sms_v0 Mobj