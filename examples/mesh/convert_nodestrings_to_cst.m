% Read nodestring of buffer area from ROSA big mesh
base_dir = '/users/modellers/rito';

% base_dir = '~/Models/FVCOM/fvcom-projects/tapas/grids/lymebay/admesh/'

project_dir = fullfile(base_dir, '/Models/FVCOM/fvcom-projects/conway/');
cd(fullfile(project_dir, 'matlab'))
addpath(fullfile(base_dir, 'Code/fvcom-toolbox/utilities'))
addpath(fullfile(base_dir, 'Code/fvcom-toolbox/fvcom_prepro/'))
%%%------------------------------------------------------------------------
%%%                          INPUT CONFIGURATION
%%%
%%% Define the model input parameters here e.g. number of tidal components,
%%% type of boundary forcing, estimated current velocity, sponge layer
%%% coefficient and radius, forcing source etc.
%%%
%%%------------------------------------------------------------------------

conf.base = fullfile(project_dir, 'run');

%%%------------------------------------------------------------------------
%%%                             Time stuff
%%%------------------------------------------------------------------------

% Case name for the model inputs and outputs
conf.casename = 'conway_v0_aqua_v16_band_nest_large_single_ND';
conf.coordinates = 'cartesian'; % 'cartesian' or 'spherical'
% conf.run_coordinates = conf.coordinates;
% Input grid UTM Zone (if applicable)



Mobj = read_sms_mesh(...
    '2dm', fullfile(conf.base, '..', 'grids',  [conf.casename, '.2dm']),...
    'coordinate', conf.coordinates, 'addCoriolis', false);

Boundary.nodes = [Mobj.read_obc_nodes{:}];
Boundary.x = Mobj.x(Boundary.nodes);
Boundary.y = Mobj.y(Boundary.nodes);

 write_SMS_cst('../grids/buffer_nodes_conway_aqua_v16.cst',Boundary.x,Boundary.y)


%% Now convert the high res subdomain nodestring to build the buffer area between coarse domain and highres domain 
 
cd(fullfile(project_dir, 'matlab'))
conf.base = fullfile(project_dir, 'run');

% Case name for the model inputs and outputs
conf.casename = 'conway_v1';

conf.coordinates = 'cartesian'; % 'cartesian' or 'spherical'

 Mobj = read_sms_mesh(...
    '2dm', fullfile(conf.base, '..', 'grids',  [conf.casename, '.2dm']),...
    'coordinate', conf.coordinates, 'addCoriolis', false);

Boundary.nodes = [Mobj.read_obc_nodes{:}];
Boundary.x = Mobj.x(Boundary.nodes);
Boundary.y = Mobj.y(Boundary.nodes);

 write_SMS_cst('../grids/boundary_nodes_conway_v1.cst',Boundary.x,Boundary.y)
