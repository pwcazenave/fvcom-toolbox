function [Plots]=do_surface_plotMatlabMap(plotOPTS,FVCOM)
%
% Function to display color maps of FVCOM variables (i.e. temperature)
%
%  [Plots]=do_surface_plot(plotOPTS,FVCOM)
%
% DESCRIPTION:
%    Generates maps of variables using m_map toolbox (patch
%
% INPUT:
%   plotOPTS   = structure array with predefined options for generating the
%   maps
%   FVCOM  = data extracted from FVCOM NC file. See read_netCDF_FVCOM for
%   details
%
%   plotOPTS.range_lat: [50.1000 50.4000]
%   plotOPTS.range_lon: [-4.5000 -3.8500]
%   plotOPTS.fig_name: 'co2_S5_slowleak'
%   plotOPTS.mesh: [1x1 struct]
%   plotOPTS.coastline_file: '../mat/tamar3_0coast.mat'
%   plotOPTS.zone: 30
%   plotOPTS.ell: 'grs80'
%   plotOPTS.var_plot: 'PH'
%   plotOPTS.clims: [6 8]
%   plotOPTS.do_mesh: 0
%   plotOPTS.nz_plot: 1
%   plotOPTS.figure: 1
%   plotOPTS.Time_record: 7.3271e+05
%   plotOPTS.nz_plot_vec: 1
%   plotOPTS.data_dec: 5
%   plotOPTS.vel_sca: 5
%   plotOPTS.pause: 0.5000

%
% OUTPUT:
%   Plots = structure array with figure handles
%
% EXAMPLE USAGE
%    [Plots]=do_surface_plot(plotOPTS,FVCOM)
%
% Author(s):
%    Ricardo Torres and Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%
%==============================================================================
%
% adds m_map to matlab paths. file is in utilities directory.
% amend according to your m_map installation paths

figure(plotOPTS.figure);
try delete(plotOPTS.PlotoutS(plotOPTS.figure).handles)
catch 
        clf
end

% generate figure with correct projection lat and lon range ellipsoid and
% zone.
if isfield(plotOPTS,'Lontick')
    MerTick = plotOPTS.Lontick;
else
    MerTick = floor(10*(diff(plotOPTS.range_lon)/5))/10;
end
if isfield(plotOPTS,'Lattick')
    ParTick = plotOPTS.Lattick;
else
    ParTick = floor(10*(diff(plotOPTS.range_lat)/5))/10;
end

axesm('mercator','MapLatLimit',plotOPTS.range_lat,'MapLonLimit',[plotOPTS.range_lon],'MeridianLabel','on',...
    'ParallelLabel','on','MLineLocation',MerTick,'PLineLocation',ParTick,'LabelUnits','dm')


% add coastline if present
if (isfield(plotOPTS,'coastline_file') && ~isempty(plotOPTS.coastline_file) )
    if isfield (plotOPTS,'PlotoutS') && ~isempty(plotOPTS.PlotoutS(plotOPTS.figure).handles)
    else
%     coast=load(plotOPTS.coastline_file);
geoshow(plotOPTS.coastline_file,'Facecolor',[0.7 .7 .7])        
    end

end
%------------------------------------------------------------------------------
% Generate maps at a give level
%------------------------------------------------------------------------------
    mstruct = gcm;


[X, Y] = mfwdtran(mstruct,plotOPTS.mesh.lat, plotOPTS.mesh.lon);

 aa=plotOPTS.Time_record;
    fprintf('Time step %i of %i\n',aa,length(FVCOM.Time_record));
    %plot map
    hold on
    ND=ndims(FVCOM.(plotOPTS.var_plot))
    switch ND
        case 2
                Plots(plotOPTS.figure).handles=patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
        'Cdata',squeeze(FVCOM.(plotOPTS.var_plot)(:,aa)),'edgecolor','interp','facecolor','interp');

        case 3
    Plots(plotOPTS.figure).handles=patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
        'Cdata',squeeze(FVCOM.(plotOPTS.var_plot)(:,plotOPTS.nz_plot,aa)),'edgecolor','interp','facecolor','interp');
    end
    caxis(plotOPTS.clims)
    colorbar h
    % check if mesh elements are required
    if plotOPTS.do_mesh
        %plot vertices
        Plots(plotOPTS.figure).mesh=patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
            'EdgeColor',[0.6 0.6 0.6],'FaceColor','none');hold on

    end
%      setm(gca,'MapLatLimit',plotOPTS.range_lat,'MapLonLimit',[plotOPTS.range_lon])
     tightmap
    %-----------------------------------------------------------------------
    % Only in my case it needs adding 6 because proj automatically determines a
    % reference long in strides of 6deg while m_map doesn't
    %------------------------------------------------------------------------
    pause(plotOPTS.pause)
%     if aa~=length(plotOPTS.Time_record)
%         delete(Plots(plotOPTS.figure).handles(:))
%     end

%  end

 
 
 
 
 
 
%  fprintf('Time step %i of %i\n',aa,length(plotOPTS.Time_record));
%     %plot map
%     hold on
%     Plots(plotOPTS.figure).handles=patch('Vertices',[lat,lon],'Faces',plotOPTS.mesh.tri,...
%         'Cdata',squeeze(FVCOM.(plotOPTS.var_plot)(:,plotOPTS.nz_plot,aa)),'edgecolor','interp','facecolor','interp');
% 
%     caxis(plotOPTS.clims)
%     colorbar
%     % check if mesh elements are required
%     if plotOPTS.do_mesh
%         %plot vertices
%         Plots(plotOPTS.figure).mesh=patch('Vertices',[lon,lat],'Faces',plotOPTS.mesh.tri,...
%             'EdgeColor',[0.6 0.6 0.6],'FaceColor','none');hold on
% 
%     end
%     setm(gca,'MapLatLimit',plotOPTS.range_lat,'MapLonLimit',[plotOPTS.range_lon])
%     %-----------------------------------------------------------------------
%     % Only in my case it needs adding 6 because proj automatically determines a
%     % reference long in strides of 6deg while m_map doesn't
%     %------------------------------------------------------------------------
%     pause(.2)
%     if aa~=length(plotOPTS.Time_record)
%         delete(Plots(plotOPTS.figure).handles(:))
%     end
% 
% end
% 
return

