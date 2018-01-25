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
clear all
matlabrc
close all
clc

base_dir = '/users/modellers/rito';
project_dir = fullfile(base_dir, '/Models/FVCOM/fvcom-projects/conway/');
cd(fullfile(project_dir, 'matlab'))

global ftbverbose
ftbverbose = 1; % be noisy

addpath(fullfile(base_dir, 'Code/fvcom-toolbox/utilities'))
addpath(fullfile(base_dir, 'Code/fvcom-toolbox/fvcom_prepro/'))
conf.obc_tides.TMD = fullfile(base_dir, 'Code/MATLAB/toolboxes/TMD/');
addpath(conf.obc_tides.TMD)
addpath(fullfile(conf.obc_tides.TMD, 'FUNCTIONS'))
addpath('scriptlets')


%%%------------------------------------------------------------------------
%%%                          INPUT CONFIGURATION
%%%
%%% Define the model input parameters here e.g. number of tidal components,
%%% type of boundary forcing, estimated current velocity, sponge layer
%%% coefficient and radius, forcing source etc.
%%%
%%%------------------------------------------------------------------------

conf.base = fullfile(project_dir, 'run');
conf.base2='';%'/users/modellers/rito/Models/git/fvcom/rosa-implementation/aqua_v12b/run'
conf.HJB_Solver_Package ='/users/modellers/rito/matlab/HJB_Solver_Package/HJB_Solver_Package/'

%%%------------------------------------------------------------------------
%%%                             Time stuff
%%%------------------------------------------------------------------------

conf.FVCOM_version = '4.0.0';

%%%------------------------------------------------------------------------
%%%                             Time stuff
%%%------------------------------------------------------------------------
% Case name for the model inputs and outputs
conf.casename = 'conway_v0';

conf.coordinates = 'cartesian'; % 'cartesian' or 'spherical'
conf.run_coordinates = 'cartesian'; % specify to what you run the model in
% conf.run_coordinates = conf.coordinates;
% Input grid UTM Zone (if applicable)
conf.utm_zone = {'30 U'}; % syntax for utm2deg


%%%------------------------------------------------------------------------
%%%                          Nesting info
%%%------------------------------------------------------------------------
%
conf.obc_type = 2;
conf.coarse.mesh = 'aqua_v16';
conf.coarse.dir = '/users/modellers/rito/Models/git/fvcom/rosa-implementation/grids/aqua_v16/';
conf.original.mesh = 'conway_v1.1';
conf.original.dir = '/users/modellers/rito/Models/FVCOM/fvcom-projects/conway/grids/';
conf.nest.mesh = 'aqua_v16_conway_v1_3_level_mesh_band';
conf.nest.dir = '/users/modellers/rito/Models/FVCOM/fvcom-projects/conway/grids/';
% Give some names to the boundaries. This must match the number of node
% strings defined in SMS. Ideally, the order of the names should match the
% order in which the boundaries were made in SMS.
conf.original.boundary_names = {'South'};
conf.coarse.boundary_names = {'West','North','East'};
conf.nest.boundary_names = {'South'};
% Nesting type (1, 2 == direct nesting, 3 == weighted)
conf.nest.type = 1;
% If we're doing type 3, we can specify the number of levels of nested
% boundaries to use. The minimum valid value is 1. For Indirect or direct
% nesting use 1
conf.nest.levels = 3; % The conway only really has 3 nested levels in common with the ROSA domain. 
                % The inner two levels had to be modified and have
                % different number of elements to the ROSA domain!! Don't
                % forget!!
conf.nest.power = 0; % zero for linear weights 1-4 increases drop in weights from initial BC position
% Ramp period for nesting inputs. Set to zero or negative to disable.



%%%------------------------------------------------------------------------
%%%                      END OF INPUT CONFIGURATION
%%%------------------------------------------------------------------------

%% Process the files needed for all months.
% Read the input mesh and bathymetry. Also creates the data necessary for
% the Coriolis correction in FVCOM.
% read original coarse mesh
outbase = fullfile(conf.base, 'input', conf.casename);
% Make the output directory if it doesn't exist.
if exist(outbase, 'dir') ~= 7
    mkdir(outbase)
end

Mobj.coarse = read_sms_mesh(...
    '2dm', fullfile(conf.coarse.dir, [conf.coarse.mesh, '.2dm']),...
    'coordinate', conf.coordinates, 'addCoriolis', false);
Mobj.nest = read_sms_mesh(...
    '2dm', fullfile(conf.nest.dir, [conf.nest.mesh, '.2dm']),...
    'coordinate', conf.coordinates, 'addCoriolis', false);

Mobj.original = read_sms_mesh(...
    '2dm', fullfile(conf.original.dir, [conf.original.mesh, '.2dm']),...
    'coordinate', conf.coordinates, 'addCoriolis', false);


% Add grid metrics.
Mobj.coarse = setup_metrics(Mobj.coarse);
Mobj.nest = setup_metrics(Mobj.nest);
Mobj.original = setup_metrics(Mobj.original);

% Parse the open boundary nodes and add accordingly.
for aa={'coarse','nest','original'}
if Mobj.(aa{1}).have_strings
    for i = 1:size(Mobj.(aa{1}).read_obc_nodes, 2)
        nodeList = double(cell2mat(Mobj.(aa{1}).read_obc_nodes(i)));
        Mobj.(aa{1}) = add_obc_nodes_list(Mobj.(aa{1}), nodeList, conf.(aa{1}).boundary_names{i}, conf.obc_type);
    end
    clear nodeList
end
end
%% Find nested range
nest.Nested_type = conf.nest.type;
nest.levels = conf.nest.levels;
nest.power = conf.nest.power;
dump = Mobj.nest;
Nested = find_nesting_region(nest, dump);

clear dump
figure(1),cla
patch('Vertices', [Mobj.original.x, Mobj.original.y], 'Faces', Mobj.original.tri, ...
    'Cdata', -Mobj.original.h, 'edgecolor', 'b', 'facecolor', 'interp');

axis('equal', 'tight')
colormap('gray')
hold on
plot(Nested.x([Nested.read_obc_nodes{:}]), Nested.y([Nested.read_obc_nodes{:}]), 'wo-')
plot(Nested.xc([Nested.read_obc_elems{:}]), Nested.yc([Nested.read_obc_elems{:}]), 'w-x')

saveas(gcf,fullfile(conf.base, ...
        'input', ...
        conf.casename,[conf.casename,'v16_nest.png']),'png')

% Plot Coarse Mesh with nesting bang 
figure(1),cla
patch('Vertices', [Nested.x, Nested.y], 'Faces', Nested.tri, ...
    'Cdata', -Nested.h, 'edgecolor', 'b', 'facecolor', 'interp');

axis('equal', 'tight')
colormap('gray')
hold on
plot(Nested.x([Nested.read_obc_nodes{:}]), Nested.y([Nested.read_obc_nodes{:}]), 'wo-')
plot(Nested.xc([Nested.read_obc_elems{:}]), Nested.yc([Nested.read_obc_elems{:}]), 'w-x')


%% Find node numbers of nested boundary nodes in coarse mesh
for nn=1:length(Nested.read_obc_nodes)
NodeList =double(Nested.read_obc_nodes{nn});
% Find the nodes that are identical in coarse domain
for ii=1:length(NodeList)
    [missmatch(ii),Mobj.coarse.nest.node{nn}(ii)]=min(abs(complex( Nested.x(NodeList(ii))-Mobj.coarse.x,Nested.y(NodeList(ii))-Mobj.coarse.y ) )); 
end
end
%% Find node numbers of nested boundary nodes in coarse mesh
for nn=1:length(Nested.read_obc_elems)
    ElemList =double(Nested.read_obc_elems{nn});
    % Find the elements that are identical in coarse domain
    for ii=1:length(ElemList)
        [missmatch(ii),Mobj.coarse.nest.elem{nn}(ii)]=min(abs(complex( Nested.xc(ElemList(ii))-Mobj.coarse.xc,Nested.yc(ElemList(ii))-Mobj.coarse.yc ) ));
    end
end


% %% Add additional nesting levels extracted from the coarse domain progressing outwards
% % I AM NOT SURE THIS IS CORRECT... IN THAT IF WE WANT MORE LEVELS WE SHOULD
% % JUST CORRECT THE HIGH RES MESH AND RE GENERATE NEST NODE LIST AND RE_RUN
% % PARENT MODEL. WITH THIS APPROACH WE WILL NEED TO GENERATE NESTING INPUTS
% % FOR HIGH RES DOMAINS WHICH WE CURRENTLY DON'T DO... ALTHOUGH IF WE WANTED
% % TO REDUCE THE NUMBER OF VERTICAL LEVELS WE WOULD HAVE TO DO ANYWAY...
% % start from external boundary 
% NodeList =double(Mobj.coarse.nest.node{1}(:));
% TR = triangulation(Mobj.coarse.tri, [Mobj.coarse.x, Mobj.coarse.y]);
% cumulative_node_idx = 1;
% cumulative_elem_idx = 1;
% for obc_idx = 1:4
%         % Given the current open boundary, find the elements connected to it
%     ti = vertexAttachments(TR, NodeList);
%     Mobj.coarse.nest.elem{end+1} = setdiff(unique([ti{:}]),[Mobj.coarse.nest.elem{:}]);
%     Mobj.coarse.nest.node{end+1} = setdiff(unique(Mobj.coarse.tri(Mobj.coarse.nest.elem{end},:)),[Mobj.coarse.nest.node{:}])';
%     NodeList=double(Mobj.coarse.nest.node{end}(:));
%     % add nest levels to conf
%     conf.nest.levels = conf.nest.levels+1;
%     % add nodes to read_obc_nodes even if they don't exist in current mesh
%     Nested.read_obc_nodes{end+1}=zeros(1,length(NodeList));
% plot(Mobj.coarse.x([Mobj.coarse.nest.node{end}]), Mobj.coarse.y([Mobj.coarse.nest.node{end}]), 'ro-')
% plot(Mobj.coarse.xc([Mobj.coarse.nest.elem{end}]), Mobj.coarse.yc([Mobj.coarse.nest.elem{end}]), 'r-x')
% end

% write nesting nodes list to be used in FVCOM coarse domain simulations
% I think the comments in mod_nesting.F are wrong and we need 
% NESTid largedomain_nodes subdomain_nodes
% for direct nesting only output boundary nodes and one interior band. This is required 
% by FVCOM as it finds the elements using bode node bands... 
% if you want relaxation nesting use more levels.


    
Nested_nodes =[Nested.read_obc_nodes{:}];
Coarse_nodes =[Mobj.coarse.nest.node{:}];
fprintf(fout,'Node_Nest Number = %d \n',length(Coarse_nodes))
for ii=1:length(Coarse_nodes)
    fprintf(fout,'%d %d %d \n',ii,Coarse_nodes(ii),Nested_nodes(ii));
end
fclose(fout)

fout=fopen(fullfile(conf.base, ...
        'input', ...
        conf.casename,[conf.casename,'latlon_v16_nest.dat']),'wt')
    
fprintf(fout,'Node_Nest Number = %d \n',length(Coarse_nodes))
for ii=1:length(Coarse_nodes)
    fprintf(fout,'%d %f %f \n',ii,Mobj.coarse.x(Coarse_nodes(ii)),Mobj.coarse.y(Coarse_nodes(ii)));
end
fclose(fout)

    
