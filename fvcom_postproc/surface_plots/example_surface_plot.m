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
%-----------------------------------------------------------
% set directories default values
%
FVCOM_root = '';
% FVCOM_data_dir='/data/medusa/rito/models/FVCOM/runCO2_leak/output/';
FVCOM_data_dir='/data/medusa/pica/models/FVCOM/runCO2_leak/output/rate_ranges/20days/';
% FVCOM_data_dir='/data/medusa/pica/models/FVCOM/runCO2_leak/output/rate_ranges/';
% FVCOM_data_dir='/data/medusa/pica/models/FVCOM/runCO2_leak/output/sponge_tests/';
FVCOM_mat_dir='../mat/';
FVCOM_plot_dir='../plots/';
FVCOM_ ='./';
%
time_offset = 678942; % from FVCOM time to matlab time
%
%------------------------------------------------------------------------------
% set casename specifics here
%------------------------------------------------------------------------------
%
[plotOPTS.range_lat ,plotOPTS.range_lon] = deal([50.1 50.4],[ -4.5 -3.85]);
base_year = datenum(2006,1,1,0,0,0);
% files_FVCOM = 'co2_S5.1.2.1_0002.nc' % Coarse, shallow, high tide, slow leak
% files_FVCOM = 'co2_S5.1.1.2_0002.nc' % Coarse, shallow, low tide, fast leak
% files_FVCOM = 'co2_V5.1.2.1_0001.nc'
% files_FVCOM = 'co2_V7.1.1.1_0001.nc'
% files_FVCOM = 'co2_V7.1.2.1_0001.nc'
% files_FVCOM = 'co2_V7.1.2.1_0002.nc'
% files_FVCOM = 'co2_S5_high_run_0001.nc';

% Currently running
% files_FVCOM = 'co2_S1_0001.nc'; % Fine, shallow, high tide, slow leak

% Scenarios
% Sponge = 0.1 i.e. low flow speeds
% files_FVCOM = 'co2_S5_high_run_fvcom_inputV5_low_flow_0001.nc';   % Coarse, fast leak
% files_FVCOM = 'co2_S5_low_run_fvcom_inputV5_low_flow_0001.nc';    % Coarse, slow leak
% files_FVCOM = 'co2_S5_pipe_run_fvcom_inputV5_low_flow_0001.nc';   % Coarse, pipe leak
% files_FVCOM = 'co2_S7_high_run_fvcom_inputV7_low_flow_0001.nc';   % Fine, fast leak
% files_FVCOM = 'co2_S7_low_run_fvcom_inputV7_low_flow_0001.nc';    % Fine, slow leak
% files_FVCOM = 'co2_S7_pipe_run_fvcom_inputV7_low_flow_0001.nc';   % Fine, pipe leak

% Sponge = 0.001 i.e. high flow speeds
% files_FVCOM = 'co2_S5_high_run_fvcom_inputV5_high_flow_0001.nc';  % Coarse, fast leak
% files_FVCOM = 'co2_S5_low_run_fvcom_inputV5_high_flow_0001.nc';   % Coarse, slow leak
% files_FVCOM = 'co2_S5_pipe_run_fvcom_inputV5_high_flow_0001.nc';  % Coarse, pipe leak
% files_FVCOM = 'co2_S7_high_run_fvcom_inputV7_high_flow.nc';       % Fine, fast leak
% files_FVCOM = 'co2_S7_low_run_fvcom_inputV7_high_flow.nc';        % Fine, slow leak
files_FVCOM = 'co2_S7_low_run_fvcom_inputV7_high_flow_0001.nc';   % Fine, slow leak
% files_FVCOM = 'co2_S7_pipe_run_fvcom_inputV7_high_flow_0001.nc';  % Fine, pipe leak

% Extras
% files_FVCOM = 'co2_S7_low_run_0001.nc'; % Fine, shallow, high tide, slow leak
% files_FVCOM = 'co2_S5_1_run_0001.nc'; % Fine, shallow, high tide, pipe day

% Sponge tests
% files_FVCOM = 'co2_S5_high_spg_1_run_fvcom_0001.nc';

if exist(fullfile(FVCOM_data_dir,files_FVCOM), 'file')~=2
    error('File not found %s',fullfile(FVCOM_data_dir,files_FVCOM))
end

casename=files_FVCOM(1:6); % 'co2_S?'
experiment_name='fastleak_hightide';
% experiment_name='slowleak_hightide';
% experiment_name='pipeleak_hightide';
% experiment_name='spongetest_0.001';


% Spring-neap
% date_range={'31/01/06 12:00:00','15/02/06 12:00:00'}; %'19/03/06 23:00:00' -1 if all available data wanted
% Neap
% date_range={'31/01/06 12:00:00','07/02/06 12:00:00'}; %'19/03/06 23:00:00' -1 if all available data wanted
% Spring
% date_range={'07/02/06 12:00:00','15/02/06 12:00:00'}; %'19/03/06 23:00:00' -1 if all available data wanted
% 4x tides at peak neaps
% date_range={'03/02/06 18:00:00','05/02/06 19:00:00'}; %'19/03/06 23:00:00' -1 if all available data wanted
% Rate testing
% date_range={'01/01/06 12:00:00','11/01/06 00:00:00'}; %'19/03/06 23:00:00' -1 if all available data wanted
% Scenarios
% date_range={'01/01/06 00:00:00','20/01/06 00:00:00'}; %'19/03/06 23:00:00' -1 if all available data wanted
% Scenarios excluding warm up
date_range={'05/01/06 23:00:00','19/01/06 00:00:00'}; %'19/03/06 23:00:00' -1 if all available data wanted

% Output figure base name
% plotOPTS.fig_name = [casename '_' experiment_name];
plotOPTS.fig_name = files_FVCOM(1:end-3);

var_2_xtractFVCOM = {'Itime','Itime2','xc','yc','h','siglay','siglev','zeta','TCO2','ph','u','v','art1','DYE'};
% var_2_xtractFVCOM = {'Itime','Itime2','xc','yc','h','siglay','siglev','zeta','TCO2','ph','art1','DYE'};
dt=1/24; % time step of output in days.
% Time record to extract
%
STnum = datenum(date_range{1},'dd/mm/yy HH:MM:SS');
ENDnum = datenum(date_range{2},'dd/mm/yy HH:MM:SS');
Time_record=STnum:dt:ENDnum;
CD=pwd;
%
%------------------------------------------------------------------------------
% Specify What indices to read from FVCOM file
%------------------------------------------------------------------------------
node_idx=-1; % to extract all nodes from netcdf file
nele_idx=-1; % to extract all elements from netcdf file
siglev_idx=-1; % to extract all water levels from netcdf file
siglay_idx=-1; % to extract all water levels from netcdf file

%
%------------------------------------------------------------------------------
% load mesh information for current casename
%------------------------------------------------------------------------------
load(fullfile(FVCOM_mat_dir, [casename, 'mesh']));
plotOPTS.mesh=mesh;

%
%------------------------------------------------------------------------------
% Specify variable and conditions of plot
%------------------------------------------------------------------------------
plotOPTS.coastline_file=[FVCOM_mat_dir 'tamar3_0coast.mat'];
plotOPTS.zone=30;
plotOPTS.var_plot = 'ph';
plotOPTS.clims=[6.0 9.0];
plotOPTS.do_mesh=0; % don't display mesh
plotOPTS.nz_plot=1;% I think is 1 for surface and 10 (the last sigma level, for bototm)
plotOPTS.Time_record=Time_record;
plotOPTS.pause=0.01; % pause between generation of plots
plotOPTS.save_output=0;

% ------------------------------------------------------------------------
%
%------------------------------------------------------------------------------
% read FVCOM data
%------------------------------------------------------------------------------
%
FVCOM_data=read_netCDF_FVCOM('time',date_range,'data_dir',FVCOM_data_dir ,...
    'file_netcdf',files_FVCOM,'node_idx',node_idx,'nele_idx',nele_idx,'siglev_idx',siglev_idx,...
    'siglay_idx',siglay_idx,'varnames',var_2_xtractFVCOM);
% put variables into a strcuture variable
for vv=1:length(var_2_xtractFVCOM)
    FVCOM.(var_2_xtractFVCOM{vv}) = FVCOM_data{vv};
end

%% Make velocity data
% eval(['vel_',files_FVCOM(1:end-3),'=sqrt(FVCOM.u(FVCOM.zeta>-4).^2+FVCOM.v(FVCOM.zeta>-4).^2);']);

%%
%------------------------------------------------------------------------------
% Do surface contour plots
%------------------------------------------------------------------------------
% select figure to plot on
plotOPTS.figure=1;
plotOPTS.Time_record=Time_record;
plotOPTS.clims=[-3 3];
PlotoutS=do_surface_plot(plotOPTS,FVCOM);
% PLotoutS has the figure handles or last handle if a time series is being
% plotted
if plotOPTS.save_output
    % Save output
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')))
    set(gcf,'renderer','painters'); % for vector output in pdfs
    print(gcf,'-dpdf','-r600',['../plots/',casename,'_',experiment_name,'_',plotOPTS.var_plot,'_contour.pdf']); % pdf
    print(gcf,'-dpng','-r600',['../plots/',casename,'_',experiment_name,'_',plotOPTS.var_plot,'_contour.png']); % png
end

%%
%------------------------------------------------------------------------------
% Do vector maps at single or multiple levels levels plots
%------------------------------------------------------------------------------
plotOPTS.figure=2;
% plotOPTS.nz_plot_vec optional. If it exists, it will overlay vectors from
% the specified layers, otherwise it uses layer from nz_plot
plotOPTS.nz_plot_vec=[1];
plotOPTS.data_dec=5;
plotOPTS.vel_sca=5; % scaling for ploting vectors with m_map m_vec.
plotOPTS.Time_record=Time_record(1:end);
% vector plots requires correct lat and lon for u and v positions: FVCOM.xc
% and FVCOM.yc. Remember to extract them
PlotoutV=do_vector_plot(plotOPTS,FVCOM)
if plotOPTS.save_output
    % Save output
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')))
    set(gcf,'renderer','painters'); % for vector output in pdfs
    print(gcf,'-dpdf','-r600',['../plots/',casename,'_',experiment_name,'_',plotOPTS.var_plot,'_sig',regexprep(num2str(plotOPTS.nz_plot_vec),'  ','_'),'.pdf']); % pdf
    print(gcf,'-dpng','-r600',['../plots/',casename,'_',experiment_name,'_',plotOPTS.var_plot,'_sig',regexprep(num2str(plotOPTS.nz_plot_vec),'  ','_'),'.png']); % png
end

%%
%------------------------------------------------------------------------------
% Do residual direction and magnitude vector maps at a single level
%------------------------------------------------------------------------------
plotOPTS.figure=3;
% plotOPTS.nz_plot_vec optional. If it exists, it will overlay vectors from
% the specified layers, otherwise it uses layer from nz_plot
plotOPTS.nz_plot_vec=[1,10];
plotOPTS.vel_sca=0.5; % scaling for ploting vectors with m_map m_vec.
% We need the model output time step (in days) for do_residual().
% Depth average the data prior residual analysis.
PlotoutR=do_residual_plot(plotOPTS,FVCOM,dt);

if plotOPTS.save_output
    % Save output
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')))
    set(gcf,'renderer','painters'); % for vector output in pdfs
    print(gcf,'-dpdf','-r600',['../plots/',casename,'_',experiment_name,'_',plotOPTS.var_plot,'_sig',regexprep(num2str(plotOPTS.nz_plot_vec),'  ','_'),'_residual.pdf']); % pdf
    print(gcf,'-dpng','-r600',['../plots/',casename,'_',experiment_name,'_',plotOPTS.var_plot,'_sig',regexprep(num2str(plotOPTS.nz_plot_vec),'  ','_'),'_residual.png']); % png
end

%% Depth averaged residual
plotOPTS.figure=4;
plotOPTS.depth_average=1; % Depth average the data prior residual analysis.
plotOPTS.vel_sca=0.5; % scaling for ploting vectors with m_map m_vec.
PlotoutR2=do_residual_plot(plotOPTS,FVCOM,dt);

if plotOPTS.save_output
    % Save output
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')))
    set(gcf,'renderer','painters'); % for vector output in pdfs
    print(gcf,'-dpdf','-r600',['../plots/',casename,'_',experiment_name,'_',plotOPTS.var_plot,'_sig0_residual.pdf']); % pdf
    print(gcf,'-dpng','-r600',['../plots/',casename,'_',experiment_name,'_',plotOPTS.var_plot,'_sig0_residual.png']); % png
end

%%
%------------------------------------------------------------------------------
% Do change in plotOPTS.var_plot (best use pH). Calculates change relative
% to base condition (i.e. first time step).
%------------------------------------------------------------------------------
plotOPTS.figure=5;
% plotOPTS.clims=[-1 0.001]; % After Jerry's ranges -- pH
plotOPTS.clims=[0 1e-5]; % DYE
plotOPTS.summed_ph=1;
plotOPTS.depth_average=1; % Depth average pH change?
plotOPTS.Time_record=Time_record(1:end);
% plotOPTS.save_intervals=[6,12,36,72,120,240];
plotOPTS.save_intervals=[];
PlotoutPH=do_ph_change_plot(plotOPTS,FVCOM,1);

if plotOPTS.save_output
    % Save output
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')))
    set(gcf,'renderer','painters'); % for vector output in pdfs
    print(gcf,'-dpdf','-r600',['../plots/',plotOPTS.fig_name,'_',plotOPTS.var_plot,'_sig0_ph_change.pdf']); % pdf
    print(gcf,'-dpng','-r600',['../plots/',plotOPTS.fig_name,'_',plotOPTS.var_plot,'_sig0_ph_change.png']); % png
end

%%
%------------------------------------------------------------------------------
% Do maximum change in pH
%------------------------------------------------------------------------------
plotOPTS.figure=6;
% plotOPTS.clims=[0 4e-3];
plotOPTS.depth_average=1; % Depth average pH change?
plotOPTS.summed_ph=0;
PlotoutPHMax=do_ph_max_change_plot(plotOPTS,FVCOM,1);

if plotOPTS.save_output
    % Save output
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')))
    set(gcf,'renderer','painters'); % for vector output in pdfs
    print(gcf,'-dpdf','-r600',['../plots/',casename,'_',experiment_name,'_',plotOPTS.var_plot,'_sig0_max_ph_change.pdf']); % pdf
    print(gcf,'-dpng','-r600',['../plots/',casename,'_',experiment_name,'_',plotOPTS.var_plot,'_sig0_max_ph_change.png']); % png
end

%%
%------------------------------------------------------------------------------
% Do volume greater or less than some threshold pH
%------------------------------------------------------------------------------
plotOPTS.figure=7;
plotOPTS.clims=[200 1.3e8];
plotOPTS.nz_plot=[1:2:10]; % 0=depth average, anything else is specified layer(s).
plotOPTS.do_mesh=0; % don't display mesh
plotOPTS.volume_thresh=+1; % -1 for less than, +1 for greater than, 0 for equal to.
[PlotoutV,totalVolume]=do_volume(plotOPTS,FVCOM,1:size(FVCOM.zeta,2),8);

%%
%------------------------------------------------------------------------------
% Do volume greater than some threshold for multiple time steps
%------------------------------------------------------------------------------
plotOPTS.threshold_change=0.2;
% .change_type: 0/1=instantaneous, any other value is hours over which to
% compare.
plotOPTS.change_type=12;
totalVolumeAffected=do_volume_change(plotOPTS,FVCOM);
fprintf('Hour %i, threshold %.2f:\t%.2fkm^{3}\n',plotOPTS.change_type,plotOPTS.threshold_change,totalVolumeAffected/1e9);


%%
%------------------------------------------------------------------------------
% Do volume greater than some threshold for multiple time steps
%------------------------------------------------------------------------------
overwrite='no';
if strcmpi(overwrite,'yes')==1

    h=1:12;
    % t=0:-0.5:-2.5;
    t=-logspace(log10(0.01),log10(2),10);
    totalVolume=nan(length(h),length(t));
    for i=h
        for j=1:length(t)
            plotOPTS.threshold_change=t(j);
            % .change_type: 0/1=instantaneous, any other value is hours over which to
            % compare.
            plotOPTS.change_type=i;
            totalVolume(i,j)=do_volume_change(plotOPTS,FVCOM);
            fprintf('Hour %i, threshold %.2f:\t%.2fkm^{3}\n',i,t(j),totalVolume(i,j)/1e9);
        end
    end

    % Save the results
%     save(['volume_analysis_h',num2str(h(1)),'-',num2str(h(end)),'_t',num2str(t(1)),'-',num2str(t(end)),'.mat'],'totalVolume','h','t');
else
    load(['volume_analysis_h',num2str(h(1)),'-',num2str(h(end)),'_t',num2str(t(1)),'-',num2str(t(end)),'.mat']);
end

%%
close all
for j=1:length(t);
    h2=h(totalVolume(:,j)~=0);
    V2=totalVolume(totalVolume(:,j)~=0);
    [slope, intercept] = logfit(h2,V2,'logy');
    yApprox = (10^intercept)*(10^slope).^h2;

    subplot(2,5,j)
    semilogy(h2,V2,'-x')
    xlabel('Hours')
    ylabel('Volume (km^{3})')
    xlim([0 12])
    ylim([1e7 max(totalVolume(:))])
    hold on
    semilogy(h2,yApprox,'kx-')
    legend('Exposed volume','Log-linear regression')
    legend('BoxOff')
    title(sprintf('Threshold %.2fpH\n',t(j)))
end


%%
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
