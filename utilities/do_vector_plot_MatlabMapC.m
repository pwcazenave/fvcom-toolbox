function [Plots]=do_vector_plot_MatlabMapC(plotOPTS,FVCOM,tt)
%
% Function to display vector maps of FVCOM currents (i.e. U,V)
%
%  [Plots]=do_vector_plot(plotOPTS,FVCOM)
%
% DESCRIPTION:
%    Generates vector maps of currents using m_map toolbox
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
%    [Plots]=do_vector_plot(plotOPTS,FVCOM)
%
% Author(s):
%    Ricardo Torres and Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%
%==============================================================================
% adds m_map to matlab paths. file is in utilities directory.
% ammend according to your m_map installation paths
figure(plotOPTS.figure);
% generate figure with correct projection lat and lon range ellipsoid and
% zone.

axesm('mercator','MapLatLimit',plotOPTS.range_lat,'MapLonLimit',[plotOPTS.range_lon],'MeridianLabel','on',...
    'ParallelLabel','on','MLineLocation',[0.1],'PLineLocation',[0.1],'LabelUnits','dm')


% add coastline if present
if (isfield(plotOPTS,'coastline_file') && ~isempty(plotOPTS.coastline_file) )
    coast=load(plotOPTS.coastline_file);
    geoshow([coast.ncst(:,2)],[coast.ncst(:,1)],'Color','black')

end
%-----------------------------------------------------------------------
% Convert element positions from FVCOM to lat and lon using m_ll2ll.m from
% utilities directory. This accesses m_map functions.
% In my case it needs adding 6 because of discrepancies between proj and m_map.
% Proj automatically determines a
% reference long in strides of 6deg while m_map doesn't
%------------------------------------------------------------------------

% only plot vectors inside lat and lon range and ...
if ~isfield(plotOPTS.mesh,'latc')
    try 
        plotOPTS.mesh.latc = nodes2elems(plotOPTS.mesh.lat,plotOPTS.mesh);
        plotOPTS.mesh.lonc = nodes2elems(plotOPTS.mesh.lon,plotOPTS.mesh);
    catch
        disp(['We have no access to nodes2elems in the fvcom_prepro directory'])
    end
end
igood = find ( plotOPTS.mesh.lonc < plotOPTS.range_lon(2) &  plotOPTS.mesh.lonc > plotOPTS.range_lon(1) &...
     plotOPTS.mesh.latc < plotOPTS.range_lat(2) &   plotOPTS.mesh.latc > plotOPTS.range_lat(1));
% decimate positions. Plot every plotOPTS.data_dec position.
igood=igood(1:plotOPTS.data_dec:end);


%------------------------------------------------------------------------
% Select how many layers to plot
%------------------------------------------------------------------------
ND=ndims(FVCOM.(plotOPTS.var_plotu))
switch ND
    case 2
        nLayers=1;
        colourSpec=[0 0 0];
        
    case 3
        if isfield(plotOPTS,'nz_plot_vec')
            nLayers=size(plotOPTS.nz_plot_vec,2);
            nLayersRange=plotOPTS.nz_plot_vec;
        else
            nLayers=size(plotOPTS.nz_plot,2);
            nLayersRange=plotOPTS.nz_plot;
        end
        % choose colors for vectors
        if nLayers==1
            colourSpec=[0 0 0];
        else
            colourSpec=colormap(hsv(nLayers));
        end
end

% Preallocate outputs
u=nan(size(igood,1),nLayers,1);
v=nan(size(igood,1),nLayers,1);

% Check if we're running
% for aa=1:length(plotOPTS.Time_record)
aa=tt;
fprintf('Time step %i of %i\n',aa,length(FVCOM.Time_record));


% Mesh goes underneath vectors.
% if plotOPTS.do_mesh
%     % plot vertices
%     Plots(plotOPTS.figure).mesh=patch('Vertices',[double(plotOPTS.mesh.lat),double(plotOPTS.mesh.lon)],'Faces',plotOPTS.mesh.tri,...
%         'EdgeColor',[0.6 0.6 0.6],'FaceColor','none');hold on
% end

switch ND
    case 2
        u(:,1)=squeeze(FVCOM.(plotOPTS.var_plotu)(igood,(aa)));
        v(:,1)=squeeze(FVCOM.(plotOPTS.var_plotv)(igood,(aa)));
    case 3
        for ii=1:nLayers
            u(:,ii)=squeeze(FVCOM.(plotOPTS.var_plotu)(igood,nLayersRange(ii),(aa)));
            v(:,ii)=squeeze(FVCOM.(plotOPTS.var_plotv)(igood,nLayersRange(ii),(aa)));
        end
end
for jj=1:nLayers
%     [Plots(plotOPTS.figure).handles{jj}]=quiver(double(plotOPTS.mesh.lonc(igood)),double(plotOPTS.mesh.latc(igood)),u,v,plotOPTS.vel_sca,'k','filled');
[Plots(plotOPTS.figure).handles{jj}]=quiverwcolorbar(double(plotOPTS.mesh.lonc(igood)),double(plotOPTS.mesh.latc(igood)),u,v,...
    plotOPTS.vel_sca,'bounds',plotOPTS.clims)
    %     [Plots(plotOPTS.figure).handles{jj}]=quiverm(double(FVCOM.latc(igood)),double(FVCOM.lonc(igood)),u./100,v./100,'k',0,'filled'),
%     [Plots(plotOPTS.figure).handles(jj),~]=m_vec(plotOPTS.vel_sca,x(igood),y(igood),...
%         u(:,jj),v(:,jj),colourSpec(jj,:),...
%         'shaftwidth',1,'headwidth',2);
    %         hold on
end
pause(plotOPTS.pause)
%     if aa~=length(plotOPTS.Time_record)
%         delete(Plots(plotOPTS.figure).handles(:))
%     end
% end
