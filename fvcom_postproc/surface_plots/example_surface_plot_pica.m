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
%    Ricardo Torres Plymouth Marine Laboratory
%
% Revision history
%   
%==============================================================================
clear
close all
clc

base = checkos;

cd([base, 'FVCOM/fvcom_co2/fvcom-toolbox/fvcom_postproc/surface_plots/'])

addpath ../../utilities
%%-----------------------------------------------------------
% set directories default values
%
FVCOM_root = '';
FVCOM_data_dir='/data/medusa/rito/models/FVCOM/runCO2_leak/output/';
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
experiment_name= 'slowleak';
[plotOPTS.range_lat ,plotOPTS.range_lon] = deal([50.1 50.4],[ -4.4 -3.85]);
base_year = datenum(2006,1,1,0,0,0);
files_FVCOM = 'co2_S5.1.2.1_0002.nc';%
date_range={'30/01/06 00:00:00','01/02/06 23:00:00'}; %'19/03/06 23:00:00' -1 if all available data wanted 
plotOPTS.fig_name = [casename '_' experiment_name];
var_2_xtractFVCOM = {'Itime','Itime2','xc','yc','h','siglay','siglev','zeta','TCO2','ph','u','v'};
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
node_idx = -1; % to extract all nodes from netcdf file
nele_idx = -1; % to extract all elements from netcdf file
siglev_idx = -1; % to extract all water levels from netcdf file
siglay_idx = -1; % to extract all water levels from netcdf file

%%
%------------------------------------------------------------------------------
% load mesh information for current casename
%------------------------------------------------------------------------------
load(fullfile(FVCOM_mat_dir, [casename, 'mesh']));
plotOPTS.mesh=mesh;

%%
%------------------------------------------------------------------------------
% Specify variable and conditions of plot
%------------------------------------------------------------------------------
plotOPTS.coastline_file=[FVCOM_mat_dir 'tamar3_0coast.mat'];
plotOPTS.zone=30;
plotOPTS.var_plot = 'ph';
plotOPTS.clims=[6.0 8.0];
plotOPTS.do_mesh=0; % 0=don't display mesh, 1=do.
plotOPTS.nz_plot=1; % I think is 1 for surface and 10 (the last sigma level, for bototm)
plotOPTS.Time_record=Time_record;
% Vector plot options
plotOPTS.data_dec=10;
plotOPTS.vel_sca=1.6;
% Optional for multiple vectors on a single plot. If missing, we'll use the
% value in plotOPTS.nz_plot.
plotOPTS.nz_plot_vec=[1,5,10]; 

%%
%------------------------------------------------------------------------------
% read FVCOM data
%------------------------------------------------------------------------------

% do extraction
FVCOM_data=read_netCDF_FVCOM_in_progress('time',date_range,'data_dir',FVCOM_data_dir ,...
    'file_netcdf',files_FVCOM,'node_idx',node_idx,'nele_idx',nele_idx,'siglev_idx',siglev_idx,...
    'siglay_idx',siglay_idx,'varnames',var_2_xtractFVCOM);
for vv=1:length(var_2_xtractFVCOM)
    FVCOM.(var_2_xtractFVCOM{vv}) = FVCOM_data{vv};
end

%% 
%------------------------------------------------------------------------------
% Do surface contour plots
%------------------------------------------------------------------------------
figure(1); clf
plotOPTS.figure=1;
Plotout=do_surface_plot(plotOPTS,FVCOM);
%%
%------------------------------------------------------------------------------
% Do vector maps at single/multiple levels plots
%------------------------------------------------------------------------------
figure(2); clf
plotOPTS.figure=2;
PlotoutV=do_vector_plot(plotOPTS,FVCOM);
% 

%%
%------------------------------------------------------------------------------
% Do residual vector plots
%------------------------------------------------------------------------------


% m_proj('UTM','lon',[ IM_range_lon],'lat',[IM_range_lat],'zon',30,'ell','grs80')
% m_grid('box','fancy')
% m_usercoast([FVCOM_mat_dir 'tamar3_0coast.mat'],'Color','k','LineWidth',3);
% %plot vertices
% %    col=0.7*[1 1 1];
% [X,Y]=m_ll2xy(mesh.lon,mesh.lat,'clip','on');
%    jnk=patch('Vertices',[X,Y],'Faces',mesh.tri,...
% 	   'EdgeColor',[0.6 0.6 0.6],'FaceColor','none');hold on
% % plot with matlab
% % HPatch=patch('Vertices',[X,Y],'Faces',mesh.tri,...
% %     'Cdata',squeeze(mesh.h),'edgecolor','interp','facecolor','interp');
% 
% for aa=1:length(Time_record)-1
% %plot map with salinity
% hold on
% HPatch=patch('Vertices',[X,Y],'Faces',mesh.tri,...
%     'Cdata',squeeze(FVCOM.(var_plot)(:,nz_plot,aa)),'edgecolor','interp','facecolor','interp');
% aa
% caxis([clims])
% colorbar
% %-----------------------------------------------------------------------
% % Only in my case it needs adding 6 because proj automatically determines a
% % reference long in strides of 6deg while m_map doesn't
% %------------------------------------------------------------------------
% pause(.2)
% delete(HPatch)
% end
% % %%
% figure(1);clf
% m_proj('UTM','lon',[ IM_range_lon],'lat',[IM_range_lat],'zon',30,'ell','grs80')
% m_grid('box','fancy')
% m_usercoast([FVCOM_mat_dir 'tamar3_0coast.mat'],'Color','k','LineWidth',3);
% [x,y]=m_ll2ll(FVCOM.xc,FVCOM.yc);x=x+6;
% 
% igood = find (x<IM_range_lon(2) & x > IM_range_lon(1) &...
%     y<IM_range_lat(2) & y>IM_range_lat(1));
% igood=igood(1:data_dec:end);
% for aa=1:length(Time_record)-1
% %plot map with salinity
% %plot vertices
% %    col=0.7*[1 1 1];
% % m_line(x+6,y,'LineStyle','none','Marker','.')
% % add vector plots on top but decimate every 200
%     u=squeeze(FVCOM.u(igood,nz_plot,(aa)));
%     v=squeeze(FVCOM.v(igood,nz_plot,(aa)));
%     [hp,ht]=m_vec(vel_sca,x(igood),y(igood),...
%         u,v,...
%         'shaftwidth',1,'headwidth',2);
%     set(hp,'facecolor','k');
%     set(hp,'edgecolor','k');
% pause(.4)
% delete(hp)
% aa
% end
% %%
% % to plot surface elevation.
% % restrict data to time range
% figure(1);clf
% m_proj('UTM','lon',[ IM_range_lon],'lat',[IM_range_lat],'zon',30,'ell','grs80')
% m_grid('box','fancy')
% m_usercoast([FVCOM_mat_dir 'tamar3_0coast.mat'],'Color','k','LineWidth',3);
% %plot vertices
% %    col=0.7*[1 1 1];
% [X,Y]=m_ll2xy(mesh.lon,mesh.lat,'clip','on');
%    jnk=patch('Vertices',[X,Y],'Faces',mesh.tri,...
% 	   'EdgeColor',[0.6 0.6 0.6],'FaceColor','none');hold on
% % plot with matlab
% % HPatch=patch('Vertices',[X,Y],'Faces',mesh.tri,...
% %     'Cdata',squeeze(mesh.h),'edgecolor','interp','facecolor','interp');
% 
% for aa=1:length(Time_record)-1
% %plot map with salinity
% hold on
% HPatch=patch('Vertices',[X,Y],'Faces',mesh.tri,...
%     'Cdata',squeeze(FVCOM.zeta(:,aa)),'edgecolor','interp','facecolor','interp');
% aa
% caxis([-2 2])
% colorbar
% %-----------------------------------------------------------------------
% % Only in my case it needs adding 6 because proj automatically determines a
% % reference long in strides of 6deg while m_map doesn't
% %------------------------------------------------------------------------
% pause(.2)
% delete(HPatch)
% end
% %%
% 
% figure(1);clf
% m_proj('UTM','lon',[ IM_range_lon],'lat',[IM_range_lat],'zon',30,'ell','grs80')
% m_grid('box','fancy')
% m_usercoast([FVCOM_mat_dir 'tamar3_0coast.mat'],'Color','k','LineWidth',3);
% [x,y]=m_ll2ll(FVCOM.xc,FVCOM.yc);x=x+6;
% 
% igood = find (x<IM_range_lon(2) & x > IM_range_lon(1) &...
%     y<IM_range_lat(2) & y>IM_range_lat(1));
% igood=igood(1:data_dec:end);
% for aa=1:length(Time_record)-1
% %plot map with salinity
% %plot vertices
% %    col=0.7*[1 1 1];
% % add vector plots on top but decimate every 200
%     u=squeeze(FVCOM.ua(igood,(aa)));
%     v=squeeze(FVCOM.va(igood,(aa)));
%     [hp,ht]=m_vec(vel_sca,x(igood),y(igood),...
%         u,v,...
%         'shaftwidth',1,'headwidth',2);
%     set(hp,'facecolor','k');
%     set(hp,'edgecolor','k');
% pause(.4)
% delete(hp)
% end
% %%
% % pcolor(squeeze(FVCOM.u(8187,:,:)));shading flat,colorbar% m_line(x+6,y,'LineStyle','none','Marker','.')
% % figure(1);clf
% % m_proj('UTM','lon',[ IM_range_lon],'lat',[IM_range_lat],'zon',30,'ell','grs80')
% % m_grid('box','fancy')
% % m_usercoast([FVCOM_mat_dir 'tamar3_0coast.mat'],'Color','k','LineWidth',3);
% % %plot vertices
% % %    col=0.7*[1 1 1];
% % [X,Y]=m_ll2xy(mesh.lon,mesh.lat,'clip','on');
% % %    jnk=patch('Vertices',[X,Y],'Faces',mesh.tri,...
% % % 	   'EdgeColor',col,'FaceColor','none');
% % % plot with matlab
% % 
% % for aa=1:length(Time_record)-1
% % %plot map with salinity
% % hold on
% % HPatch=patch('Vertices',[X,Y],'Faces',mesh.tri,...
% %     'Cdata',squeeze(FVCOM.(var_plot)(:,nz_plot,aa)),'edgecolor','interp','facecolor','interp');
% % aa
% % caxis([clims])
% % colorbar
% % %-----------------------------------------------------------------------
% % % Only in my case it needs adding 6 because proj automatically determines a
% % % reference long in strides of 6deg while m_map doesn't
% % %------------------------------------------------------------------------
% % pause(0.4)
% % delete(HPatch)
% % end
% % figure(1);clf
% % m_proj('UTM','lon',[ IM_range_lon],'lat',[IM_range_lat],'zon',30,'ell','grs80')
% % m_grid('box','fancy')
% % m_usercoast([FVCOM_mat_dir 'tamar3_0coast.mat'],'Color','k','LineWidth',3);
% % [x,y]=m_ll2ll(FVCOM.xc,FVCOM.yc);x=x+6;
% % 
% % igood = find (x<IM_range_lon(2) & x > IM_range_lon(1) &...
% %     y<IM_range_lat(2) & y>IM_range_lat(1));
% % igood=igood(1:data_dec:end);
% % for aa=1:length(Time_record)-1
% % %plot map with salinity
% % %plot vertices
% % %    col=0.7*[1 1 1];
% % % m_line(x+6,y,'LineStyle','none','Marker','.')
% % % add vector plots on top but decimate every 200
% %     u=squeeze(FVCOM.u(igood,nz_plot,(aa)));
% %     v=squeeze(FVCOM.v(igood,nz_plot,(aa)));
% %     [hp,ht]=m_vec(vel_sca,x(igood),y(igood),...
% %         u,v,...
% %         'shaftwidth',1,'headwidth',2);
% %     set(hp,'facecolor','k');
% %     set(hp,'edgecolor','k');
% % pause(0.5)
% % delete(hp)
% % end