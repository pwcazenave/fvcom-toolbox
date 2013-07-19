% Sample script to extract and generate transect plots of tracer variables 
% Modify to suit your requirements
%
% DESCRIPTION:
%    Extracts data from FVCOM and generates vertical sliced contour plots 
%
% INPUT: 
%   
%   
%           
% OUTPUT:
%   
% EXAMPLE USAGE
%    See scripts and modify at will 
%
% Author(s):  
%    Ricardo Torres Plymouth Marine Laboratory
%
% Revision history
%   
%==============================================================================

clear
close all
clc

addpath ../../utilities
%%-----------------------------------------------------------
% set directories default values
%
FVCOM_root = '';
FVCOM_data_dir='/home_nfs/rito/models/FVCOM/runCO2_leak/output/';
FVCOM_mat_dir='../mat/';
FVCOM_plot_dir='../plots/';
FVCOM_ ='./';
%
time_offset = 678942; % from FVCOM time to matlab time
%%
%------------------------------------------------------------------------------
% set casename specifics here
%------------------------------------------------------------------------------
%
casename='co2_S5';
experiment_name='slowleak';
base_year = datenum(2006,1,1,0,0,0);
files_FVCOM ='co2_S5.1.2.1_0002.nc';%
date_range={'30/01/06 00:00:00','01/02/06 23:00:00'}; %'19/03/06 23:00:00' -1 if all available data wanted 
plotOPTS.fig_name = [casename '_' experiment_name];
var_2_xtractFVCOM = {'Itime','Itime2','xc','yc','h','siglay','siglev','zeta','TCO2','ph','u','v','temp','density'};
dt=1/24; % time step of output in days.
% Time record to extract
%
STnum = datenum(date_range{1},'dd/mm/yy HH:MM:SS');
ENDnum = datenum(date_range{2},'dd/mm/yy HH:MM:SS');
Time_record=STnum:dt:ENDnum;
CD=pwd;
%%
%------------------------------------------------------------------------------
% Specify What indices to read from FVCOM file
%------------------------------------------------------------------------------
siglev_idx=-1; % to extract all water levels from netcdf file
siglay_idx=-1; % to extract all water levels from netcdf file

%------------------------------------------------------------------------------
load(fullfile(FVCOM_mat_dir, [casename, 'mesh']));
plotOPTS.mesh=mesh;
%------------------------------------------------------------------------------
% select transect nodes with GUI
%------------------------------------------------------------------------------
trn_nodes= transect_nodes_screen(mesh,'nodes')
node_idx=trn_nodes.idx;
nele_idx=-1;
%------------------------------------------------------------------------------
% calculate distance along transect
% ------------------------------------------------------------------------
trn_dist=[trn_nodes.x(1)-trn_nodes.x(:),...
    trn_nodes.y(1)-trn_nodes.y(:)];
trn_dis=cumsum(abs(complex(trn_dist(:,1),trn_dist(:,2))));
%
%%
%------------------------------------------------------------------------------
% Specify variable and conditions of plot
%------------------------------------------------------------------------------
plotOPTS.var_plot = 'ph';
plotOPTS.clims=[0.026418 0.026419];
plotOPTS.do_mesh=0; % don't display mesh
plotOPTS.trn_dis=trn_dis;
plotOPTS.trn_nodes=trn_nodes;
[plotOPTS.range_lat ,plotOPTS.range_lon] = deal([50.1 50.4],[ -4.5 -3.85]);
plotOPTS.coastline_file=[FVCOM_mat_dir 'tamar3_0coast.mat'];
plotOPTS.zone=30;
plotOPTS.ell='grs80';
plotOPTS.do_mesh=0; % don't display mesh

%% ------------------------------------------------------------------------
% 
%------------------------------------------------------------------------------
% read FVCOM data
%------------------------------------------------------------------------------
% 
FVCOM_data=read_netCDF_FVCOM('time',date_range,'data_dir',FVCOM_data_dir ,...
    'file_netcdf',files_FVCOM,'node_idx',node_idx,'nele_idx',nele_idx,'siglev_idx',siglev_idx,...
    'siglay_idx',siglay_idx,'varnames',var_2_xtractFVCOM)
% put variables into a strcuture variable
for vv=1:length(var_2_xtractFVCOM)
    FVCOM.(var_2_xtractFVCOM{vv}) = FVCOM_data{vv};
end

cd(CD)
plotOPTS.figure=1;
[Plots]=do_transect_plot(plotOPTS,FVCOM)



