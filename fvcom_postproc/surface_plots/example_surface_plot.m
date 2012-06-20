% Sample script to extract and generate m_map contours of tracer variables 
% Modify to suit your requirements
%
% DESCRIPTION:
%    Extracts data from FVCOM and generates contour plots using m_map 
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
%    Ricardo Torres and Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   
%==============================================================================
addpath ../utilities
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
experiment_name='slowleak'
[plotOPTS.range_lat ,plotOPTS.range_lon] = deal([50.1 50.4],[ -4.5 -3.85]);
base_year = datenum(2006,1,1,0,0,0);
files_FVCOM ='co2_S5.1.2.1_0002.nc'%
date_range={'30/01/06 00:00:00','01/02/06 23:00:00'}; % -1 if all available data wanted 
plotOPTS.fig_name = [casename '_' experiment_name];
var_2_xtractFVCOM = {'Itime','Itime2','xc','yc','h','siglay','siglev','zeta','salinity','u','v'};
dt=1/24; % time step of output in days.
% Time record to extract
%
STnum = datenum(date_range{1},'dd/mm/yy HH:MM:SS');
ENDnum = datenum(date_range{2},'dd/mm/yy HH:MM:SS');
Time_record=STnum:dt:ENDnum;
CD=pwd
%%
%------------------------------------------------------------------------------
% Specify What indices to read from FVCOM file
%------------------------------------------------------------------------------
node_idx= -1 % to extract all nodes from netcdf file
nele_idx=-1 % to extract all elements from netcdf file
siglev_idx=-1 % to extract all water levels from netcdf file
siglay_idx=-1 % to extract all water levels from netcdf file

%%
%------------------------------------------------------------------------------
% load mesh information for current casename This is generated at preproc
% stage (see %%%% matlab script)
%------------------------------------------------------------------------------
load(fullfile(FVCOM_mat_dir, [casename, 'mesh']));
plotOPTS.mesh=mesh;

%%
%------------------------------------------------------------------------------
% Specify variable and conditions of plot
%------------------------------------------------------------------------------
% if you have a coastline that can be directly used by m_map here is the
% place to put it!
plotOPTS.coastline_file=[FVCOM_mat_dir 'tamar3_0coast.mat'];
% m_map projection information
plotOPTS.zone=30;
plotOPTS.ell='grs80';
% variable to plot
plotOPTS.var_plot = 'PH';
% color limits for the colorbar
plotOPTS.clims=[6.0 8.0]
plotOPTS.do_mesh=0; % don't display mesh 1 to overlay mesh
plotOPTS.nz_plot=1;% layer to plot
plotOPTS.Time_record=Time_record; % time steps to make plots (1 for single or more for animation)
%% ------------------------------------------------------------------------
%% 
%------------------------------------------------------------------------------
% read FVCOM data. See read_netCDF_FVCOM.m for help 
%------------------------------------------------------------------------------
% 
FVCOM_data=read_netCDF_FVCOM('time',date_range,'data_dir',FVCOM_data_dir ,...
    'file_netcdf',files_FVCOM,'node_idx',node_idx,'nele_idx',nele_idx,'siglev_idx',siglev_idx,...
    'siglay_idx',siglay_idx,'varnames',var_2_xtractFVCOM)
% put variables into a strcuture variable
for vv=1:length(var_2_xtractFVCOM)
    FVCOM.(var_2_xtractFVCOM{vv}) = FVCOM_data{vv};
end

%% 
%------------------------------------------------------------------------------
% Do surface contour plots
%------------------------------------------------------------------------------
% select figure to plot on
plotOPTS.figure=1;
PLotoutS=do_surface_plot(plotOPTS,FVCOM)
% PLotoutS has the figure handles or last handle if a time series is being
% plot
%------------------------------------------------------------------------------
% Do vector maps at single or multiple levels plots
%------------------------------------------------------------------------------
plotOPTS.figure=2;
% plotOPTS.nz_plot_vec optional. If it exists, it will overlay vectors from
% the specified layers, otherwise it uses layer from nz_plot
plotOPTS.nz_plot_vec=[1,10];
plotOPTS.data_dec=5;
plotOPTS.vel_sca=5; % scaling for ploting vectors with m_map m_vec.
plotOPTS.pause=0.5 % pause between generation of plots
% vector plots requires correct lat and lon for u and v positions: FVCOM.xc
% and FVCOM.yc. Remember to extract them
PLotoutV=do_vector_plot(plotOPTS,FVCOM)
