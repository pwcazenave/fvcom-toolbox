% Get transects for the coarse and fine grids for the pH profiles

% Load the results and the mesh

clear
close all
clc

% Redo the transect or load a saved one?
reselecttransect='no';

addpath ../../utilities
addpath /users/modellers/rito/Models/fvcom-toolbox/utilities/
%-----------------------------------------------------------
% set directories default values
%
FVCOM_root = '';
% FVCOM_data_dir='/data/medusa/rito/models/FVCOM/runCO2_leak/output/';
FVCOM_data_dir='/data/medusa/pica/models/FVCOM/runCO2_leak/output/scenarios/';
% FVCOM_data_dir='/tmp/data/20days/';
% FVCOM_data_dir='/data/medusa/pica/models/FVCOM/runCO2_leak/output/rate_ranges/';
% FVCOM_data_dir='/data/medusa/pica/models/FVCOM/runCO2_leak/output/sponge_tests/';
FVCOM_mat_dir='/tmp/fvcom_co2/mat/';
plotOPTS.FVCOM_plot_dir='/tmp/fvcom_co2/plots/';
FVCOM_ ='./';
%
time_offset = 678942; % from FVCOM time to matlab time
%
%------------------------------------------------------------------------------
% set casename specifics here
%------------------------------------------------------------------------------
[plotOPTS.range_lat ,plotOPTS.range_lon] = deal([50.2 50.4],[ -4.34 -3.91]);
base_year = datenum(2006,1,1,0,0,0);

% Use two results only -- doesn't matter what co2 scenario, we're only
% interested in nodes and z
allFiles={
    'co2_S5_high_run_fvcom_inputV5_low_flow_0001.nc',...   % Coarse, fast leak
    'co2_S7_low_run_fvcom_inputV7_low_flow_0001.nc'...    % Fine, slow leak
};

for fileIdx=1%:size(allFiles,2)
    
    files_FVCOM = allFiles{fileIdx};

    if ~isempty(strfind(files_FVCOM,'S5')) || ~isempty(strfind(files_FVCOM,'V5'))
        saveFile=fullfile(FVCOM_mat_dir,'S5_east-west_transect.mat');
    elseif ~isempty(strfind(files_FVCOM,'S7')) || ~isempty(strfind(files_FVCOM,'V7'))
        saveFile=fullfile(FVCOM_mat_dir,'S7_east-west_transect.mat');
    else
        error('Unknown grid resolution type. Can''t set output file name.')
    end

    
    if exist(fullfile(FVCOM_data_dir,files_FVCOM), 'file')~=2
        error('File not found %s',fullfile(FVCOM_data_dir,files_FVCOM))
    end
    
    % Scenarios excluding warm up
    date_range={'05/01/06 22:00:00','19/01/06 00:00:00'}; %'19/03/06 23:00:00' -1 if all available data wanted 

    casename=files_FVCOM(1:6); % 'co2_S?'

    % Output figure base name
    % plotOPTS.fig_name = [casename '_' experiment_name];
    plotOPTS.fig_name = files_FVCOM(1:end-3);

%     var_2_xtractFVCOM = {'Itime','Itime2','xc','yc','h','siglay','siglev','zeta','TCO2','ph','art1','DYE'};
    var_2_xtractFVCOM = {'Itime','Itime2','xc','yc','h','siglay','siglev','zeta','TCO2','ph','art1','DYE','u','v'};
    dt=1/24; % time step of output in days.

    % Time record to extract
    %
    STnum = datenum(date_range{1},'dd/mm/yy HH:MM:SS');
    ENDnum = datenum(date_range{2},'dd/mm/yy HH:MM:SS');
    Time_record=STnum:dt:ENDnum;
    CD=pwd;
    %
    %------------------------------------------------------------------------------
    % Specify what indices to read from FVCOM file
    %------------------------------------------------------------------------------
    siglev_idx=-1; % to extract all water levels from netcdf file
    siglay_idx=-1; % to extract all water levels from netcdf file

    %
    %------------------------------------------------------------------------------
    % load mesh information for current casename
    %------------------------------------------------------------------------------
    load(fullfile(FVCOM_mat_dir, [casename, 'mesh']));
    plotOPTS.mesh=mesh;

    %% Select the transect
    if strcmpi(reselecttransect,'yes')
        trn_nodes=transect_nodes_screen(mesh,'nodes')
    else
        trn_nodes=load(saveFile);
        trn_nodes=trn_nodes.trn_nodes;
    end
    node_idx=trn_nodes.idx;
    nele_idx=-1;
    
    %
    %------------------------------------------------------------------------------
    % Specify variable and conditions of plot
    %------------------------------------------------------------------------------
    plotOPTS.coastline_file=[FVCOM_mat_dir 'tamar3_0coast.mat'];
    plotOPTS.zone=30;
    plotOPTS.Time_record=Time_record;
    plotOPTS.ell='grs80';
    plotOPTS.do_mesh=0; % don't display mesh
    plotOPTS.var_plot = 'ph';
    plotOPTS.trn_nodes=trn_nodes;

    % Get distance along line
    xNorm=trn_nodes.x-min(trn_nodes.x);
    yNorm=trn_nodes.y-min(trn_nodes.y);
    xDist=sqrt(xNorm.^2+yNorm.^2);
    plotOPTS.trn_dis=xDist;    
    
    %------------------------------------------------------------------------------
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


    
    % Check the position of the transect
%     cd(CD)
%     plotOPTS.figure=1;
%     plotOPTS.trn_nodes=trn_nodes;
%     [Plots]=do_transect_plot(plotOPTS,FVCOM);

    % Check the depths we've extracted
    figure
    plot(xDist,-FVCOM.h)
    
    % Save the transect information to the save file
    if strcmpi(reselecttransect,'yes')
        save(saveFile,'trn_nodes')
    end
    
    %% Save a figure of the transect location
    
    % Re-read the results
    %------------------------------------------------------------------------------
    % 
    %------------------------------------------------------------------------------
    % read FVCOM data
    %------------------------------------------------------------------------------
    % 
    node_idx_all=-1; % to extract all nodes from netcdf file
    nele_idx_all=-1; % to extract all elements from netcdf file

    FVCOM_data=read_netCDF_FVCOM('time',date_range,'data_dir',FVCOM_data_dir ,...
        'file_netcdf',files_FVCOM,'node_idx',node_idx_all,'nele_idx',nele_idx_all,'siglev_idx',siglev_idx,...
        'siglay_idx',siglay_idx,'varnames',var_2_xtractFVCOM);
    % put variables into a strcuture variable
    for vv=1:length(var_2_xtractFVCOM)
        FVCOM.(var_2_xtractFVCOM{vv}) = FVCOM_data{vv};
    end

    plotOPTS.clims=[-30 0];
    m_mappath;
    figure(100); clf
    m_proj('UTM','lon',[plotOPTS.range_lon],'lat',[plotOPTS.range_lat],'zon',plotOPTS.zone,'ell','grs80')
    m_grid('box','fancy')
    m_usercoast(plotOPTS.coastline_file,'Color','k','LineWidth',3);
    [X,Y]=m_ll2xy(plotOPTS.mesh.lon,plotOPTS.mesh.lat,'clip','on');
    hold on
    patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
        'Cdata',-FVCOM.h,'edgecolor','interp','facecolor','interp');
    caxis(plotOPTS.clims)
    ch=colorbar;
    colormap(jet)
    ylabel(ch,'Depth (m)')
    plot(X(node_idx), Y(node_idx),'-w.','LineWidth',2,'MarkerSize',5)
%     plot(sort(X(node_idx)), sort(Y(node_idx),1,'descend'),'-w.','LineWidth',2,'MarkerSize',5)
    
    fprintf('Saving transect location figure... ')
    set(findobj(gcf,'Type','text'),'FontSize',10)
    %set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize'))) 
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'renderer','painters'); % for vector output in pdfs
    print(gcf,'-dpdf','-r600',[plotOPTS.FVCOM_plot_dir,plotOPTS.var_plot,'/pdf/',casename,'_transect_location.pdf']); % pdf
    %print(gcf,'-dpng','-r600',[plotOPTS.FVCOM_plot_dir,plotOPTS.var_plot,'/png/',plotOPTS.fig_name,'_layer=',num2str(plotOPTS.nz_plot),'_',plotOPTS.var_plot,'_change_',num2str(plotOPTS.save_intervals(aa)),'.png']); % png
    fprintf('done.\n')

    
    
end

