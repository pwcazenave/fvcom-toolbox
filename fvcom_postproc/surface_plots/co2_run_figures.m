% Script to make Jerry's figures. Hopefully with as little commenting in
% and out as possible.
% 
% DESCRIPTION:
%   Extracts the data from the FVCOM output and generates a range of plots.
% 
% INPUT: 
%   FVCOM model output files (*.nc)
% 
% OUTPUT:
%   Figures in FVCOM_plot_dir
% 
% Author(s):
%   Pierre Cazenave Plymouth Marine Laboratory
% 
% Revision history
%   2012/05/28 - Created.
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
FVCOM_data_dir='/data/medusa/pica/models/FVCOM/runCO2_leak/output/scenarios/';
% FVCOM_data_dir='/tmp/data/20days/';
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
[plotOPTS.range_lat ,plotOPTS.range_lon] = deal([50.2 50.4],[ -4.34 -3.91]);
plotOPTS.FVCOM_plot_dir=FVCOM_plot_dir;
base_year = datenum(2006,1,1,0,0,0);

fout=fopen([FVCOM_plot_dir, 'volumes_temp.csv'],'w');
 % add header
fprintf(fout,'%s\n','filename,"volume (m3) -0.2","volume (m3) -0.5","volume (m3) -1.0"');

allFiles={
    % Coarse
    % Low -> high -> med flow
    % High -> low -> pipe
    'co2_S5_high_run_fvcom_inputV5_low_flow_0001.nc',...   % Coarse, fast leak
    'co2_S5_low_run_fvcom_inputV5_low_flow_0001.nc',...    % Coarse, slow leak
    'co2_S5_pipe_run_fvcom_inputV5_low_flow_0001.nc',...   % Coarse, pipe leak
    'co2_S5_high_run_fvcom_inputV5_high_flow_0001.nc',...  % Coarse, fast leak
    'co2_S5_low_run_fvcom_inputV5_high_flow_0001.nc',...   % Coarse, slow leak
    'co2_S5_pipe_run_fvcom_inputV5_high_flow_0001.nc',...  % Coarse, pipe leak
    'co2_S5_high_run_fvcom_inputV5_mid_flow_0001.nc',...  % Coarse, fast leak
    'co2_S5_low_run_fvcom_inputV5_mid_flow_0001.nc',...   % Coarse, slow leak
    'co2_S5_pipe_run_fvcom_inputV5_mid_flow_0001.nc',...  % Coarse, pipe leak
    ... % Fine
    ... % Low -> high -> med flow. 
    ... % High -> low -> pipe
    ...%'co2_S7_high_run_fvcom_inputV7_low_flow_0001.nc',...   % Fine, fast leak
    'co2_S7_low_run_fvcom_inputV7_low_flow_0001.nc',...    % Fine, slow leak
    ...%'co2_S7_pipe_run_fvcom_inputV7_low_flow_0001.nc',...   % Fine, pipe leak
    'co2_S7_high_run_fvcom_inputV7_high_flow_0001.nc',...  % Fine, fast leak
    ...%'co2_S7_low_run_fvcom_inputV7_high_flow_0001.nc',...   % Fine, slow leak
    ...%'co2_S7_pipe_run_fvcom_inputV7_high_flow_0001.nc',...  % Fine, pipe leak
    ...%'co2_S7_high_run_fvcom_inputV7_mid_flow_0001.nc',...  % Fine, fast leak
    'co2_S7_low_run_fvcom_inputV7_mid_flow_0001.nc'...   % Fine, slow leak
    ...%'co2_S7_pipe_run_fvcom_inputV7_mid_flow_0001.nc'  % Fine, pipe leak
};

% Scenarios
% Sponge = 0.1 i.e. slow flow speeds
% allFiles = {'co2_S5_high_run_fvcom_inputV5_low_flow_0001.nc'};   % Coarse, fast leak
% allFiles = {'co2_S5_low_run_fvcom_inputV5_low_flow_0001.nc'};    % Coarse, slow leak
% allFiles = {'co2_S5_pipe_run_fvcom_inputV5_low_flow_0001.nc'};   % Coarse, pipe leak
% allFiles = {'co2_S7_high_run_fvcom_inputV7_low_flow_0001.nc'};   % Fine, fast leak
% allFiles = {'co2_S7_low_run_fvcom_inputV7_low_flow_0001.nc'};    % Fine, slow leak
% allFiles = {'co2_S7_pipe_run_fvcom_inputV7_low_flow_0001.nc'};   % Fine, pipe leak

% Sponge = 0.001 i.e. fast flow speeds
% allFiles = {'co2_S5_high_run_fvcom_inputV5_high_flow_0001.nc'};  % Coarse, fast leak
% allFiles = {'co2_S5_low_run_fvcom_inputV5_high_flow_0001.nc'};   % Coarse, slow leak
% allFiles = {'co2_S5_pipe_run_fvcom_inputV5_high_flow_0001.nc'};  % Coarse, pipe leak
% allFiles = {'co2_S7_high_run_fvcom_inputV7_high_flow_0001.nc'};  % Fine, fast leak
% allFiles = {'co2_S7_low_run_fvcom_inputV7_high_flow_0001.nc'};   % Fine, slow leak
% allFiles = {'co2_S7_pipe_run_fvcom_inputV7_high_flow_0001.nc'};  % Fine, pipe leak

% Sponge = 0.01 i.e. medium flow speeds
% allFiles = {'co2_S5_high_run_fvcom_inputV5_mid_flow_0001.nc'};  % Coarse, fast leak
% allFiles = {'co2_S5_low_run_fvcom_inputV5_mid_flow_0001.nc'};   % Coarse, slow leak
% allFiles = {'co2_S5_pipe_run_fvcom_inputV5_mid_flow_0001.nc'};  % Coarse, pipe leak
% allFiles = {'co2_S7_high_run_fvcom_inputV7_mid_flow_0001.nc'};  % Fine, fast leak
% allFiles = {'co2_S7_low_run_fvcom_inputV7_mid_flow_0001.nc'};   % Fine, slow leak
% allFiles = {'co2_S7_pipe_run_fvcom_inputV7_mid_flow_0001.nc'};  % Fine, pipe leak

for fileIdx=1:size(allFiles,2)
    
    files_FVCOM = allFiles{fileIdx};

    if ~isempty(strfind(files_FVCOM,'high_run'))
        plotOPTS.clims=[-0.001 0]; % High leak
        plotOPTS.vlims=[8.16 8.165];
    elseif ~isempty(strfind(files_FVCOM,'low_run'))
        plotOPTS.clims=[-0.000005 0]; % Low leak
        plotOPTS.vlims=[8.16248 8.16253];
    elseif ~isempty(strfind(files_FVCOM,'pipe_run'))
        plotOPTS.clims=[-2 0]; % Pipe leak
        plotOPTS.vlims=[6 8.5];
    else
        error('Unknown model type. Can''t set colour range.')
    end
    
    if ~isempty(strfind(files_FVCOM,'S5')) || ~isempty(strfind(files_FVCOM,'V5'))
        transectPoints=load(fullfile(FVCOM_mat_dir,'S5_east-west_transect.mat'));
    elseif ~isempty(strfind(files_FVCOM,'S7')) || ~isempty(strfind(files_FVCOM,'V7'))
        transectPoints=load(fullfile(FVCOM_mat_dir,'S5_east-west_transect.mat'));
    else
        error('Unknown grid resolution type. Can''t import saved transect.')
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
    plotOPTS.Time_record=Time_record;
    plotOPTS.pause=0.01; % pause between generation of plots

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

    %------------------------------------------------------------------------------
    %
    %------------------------------------------------------------------------------
    % get the velocity stats
    %------------------------------------------------------------------------------
    %
    [velMin,velMax,velMed,velMean,velStd]=get_velocity(FVCOM.u,FVCOM.v);
    fprintf('\n\n%s,%4f,%4f,%4f,%4f,%4f\n\n',files_FVCOM,velMin,velMax,velMed,velMean,velStd)


    %%
    %------------------------------------------------------------------------------
    % Do change in plotOPTS.var_plot (e.g. pH). Calculates change relative
    % to base condition (i.e. first time step).
    %------------------------------------------------------------------------------
    if 0
        plotOPTS.figure=1;
        plotOPTS.var_plot='ph';
        plotOPTS.do_mesh=0; % Add the mesh?

        % not doing this:
        % plotOPTS.summed_ph=1; % Cumulative (1) or instantaneous (0) 

        plotOPTS.depth_average=0; % Depth average pH change?
        plotOPTS.nz_plot_range=[1,10]; % If not depth averaged, which sigma layer(s)?

        plotOPTS.save_intervals=[6,12,36,72,120,240]; % Specified intervals only
        plotOPTS.Time_record=Time_record(plotOPTS.save_intervals);
        % plotOPTS.save_intervals=[];
        % plotOPTS.Time_record=Time_record;

        if ~isempty(strfind(files_FVCOM,'low_run')) || ~isempty(strfind(files_FVCOM,'high_run'))
            % We're doing low/high flow, so set the alternate colour
            % palette flag
            plotOPTS.altColours=1;
        else
            plotOPTS.altColours=0;
        end
        
        plotOPTS.save_output=1;

        for ii=1:size(plotOPTS.nz_plot_range,2)
            plotOPTS.nz_plot=plotOPTS.nz_plot_range(ii);
            PlotoutPH=do_ph_change_plot(plotOPTS,FVCOM,1);
        end
    end
    %%
    %------------------------------------------------------------------------------
    % Do maximum change in plotOPTS.var_plot (e.g. pH) across all time steps.
    %------------------------------------------------------------------------------
    if 0
        plotOPTS.figure=2;
        plotOPTS.var_plot='ph';
        plotOPTS.do_mesh=0; % Add the mesh?

        % not doing this:
        plotOPTS.summed_ph=0; % Cumulative (1) or instantaneous (0)
        
        plotOPTS.depth_average=0; % Depth average pH change?
        plotOPTS.nz_plot_range=[1,10]; % If not depth averaged, which sigma layer(s)?
        plotOPTS.Time_record=Time_record;

        if ~isempty(strfind(files_FVCOM,'low_run')) || ~isempty(strfind(files_FVCOM,'high_run'))
            % We're doing low/high flow, so set the alternate colour
            % palette flag
            plotOPTS.altColours=1;
        else
            plotOPTS.altColours=0;
        end
        
        plotOPTS.save_output=1;

        for ii=1:size(plotOPTS.nz_plot_range,2)
            plotOPTS.nz_plot=plotOPTS.nz_plot_range(ii);
            PlotoutPHMax=do_ph_max_change_plot(plotOPTS,FVCOM,1);
        end    
    end
    %%
    %------------------------------------------------------------------------------
    % Do vertical profile along a transect
    %------------------------------------------------------------------------------
    if 1
        plotOPTS.figure=3;
        plotOPTS.var_plot='ph';
        plotOPTS.do_mesh=0; % Add the mesh?
        plotOPTS.save_intervals=[6,12,36,72,120,240]; % Specified intervals only
        plotOPTS.Time_record=Time_record(plotOPTS.save_intervals);

        if ~isempty(strfind(files_FVCOM,'low_run')) || ~isempty(strfind(files_FVCOM,'high_run'))
            % We're doing low/high flow, so set the alternate colour
            % palette flag
            plotOPTS.altColours=1;
        else
            % Otherwise we'll just use the jet colour scheme
            plotOPTS.altColours=0;
        end
        
        plotOPTS.save_output=1;
    
        PlotoutPH=do_ph_vertical_profile(plotOPTS,FVCOM,transectPoints);
    end
    %%
    %------------------------------------------------------------------------------
    % Do volume greater than some threshold for multiple time steps
    %------------------------------------------------------------------------------
    if 0
        plotOPTS.threshold_change_range=[-0.2,-0.5,-1.0];
        plotOPTS.var_plot='ph';
        plotOPTS.nz_plot=1; % Do seabed only
        % .change_type: 0=instantaneous, any other value is hours over which to
        % compare.
        plotOPTS.change_type=0;
        fprintf(fout,'%s,',files_FVCOM);
        for ii=1:size(plotOPTS.threshold_change_range,2)
            plotOPTS.threshold_change=plotOPTS.threshold_change_range(ii);
            totalVolumeAffected=do_volume_change(plotOPTS,FVCOM);
            % For the second version, give a background condition (i.e. pre-leak).
            % Also, this won't do the n-hour check for continuous runs.
            % do_volume_change does, but does so wrong (sums volume incorrectly, I
            % think).
            totalVolumeAffected2=do_volume_change2(plotOPTS,FVCOM,1);
        %     fprintf('Hour %i, threshold %.2f:\t%.2f km^{3}\n',plotOPTS.change_type,plotOPTS.threshold_change,totalVolumeAffected/1e9);
        %     fprintf('%.2f\n',totalVolumeAffected/1e9);
            fprintf(fout,'%.2f,',totalVolumeAffected2); % output in m3
        end
        fprintf(fout,'\n');
    end
end
fclose(fout);