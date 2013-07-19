function [Plots]=do_vector_plot(plotOPTS,FVCOM)
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
m_mappath;
% adds m_map to matlab paths. file is in utilities directory.
% ammend according to your m_map installation paths
figure(plotOPTS.figure);
% generate figure with correct projection lat and lon range ellipsoid and
% zone.
m_proj('UTM','lon',[plotOPTS.range_lon],'lat',[plotOPTS.range_lat],'zon',plotOPTS.zone,'ell',plotOPTS.ell)
m_grid('box','fancy')
% add coastline if present
if (isfield(plotOPTS,'coastline_file') && ~isempty(plotOPTS.coastline_file) )
    m_usercoast(plotOPTS.coastline_file,'Color','k','LineWidth',3);
end
%-----------------------------------------------------------------------
% Convert element positions from FVCOM to lat and lon using m_ll2ll.m from
% utilities directory. This accesses m_map functions.
% In my case it needs adding 6 because of discrepancies between proj and m_map.
% Proj automatically determines a
% reference long in strides of 6deg while m_map doesn't
%------------------------------------------------------------------------

[x,y]=m_ll2ll(FVCOM.xc,FVCOM.yc); x=x+6;

% only plot vectors inside lat and lon range and ...
igood = find (x < plotOPTS.range_lon(2) & x > plotOPTS.range_lon(1) &...
    y < plotOPTS.range_lat(2) & y > plotOPTS.range_lat(1));
% decimate positions. Plot every plotOPTS.data_dec position.
igood=igood(1:plotOPTS.data_dec:end);
%------------------------------------------------------------------------
% Select how many layers to plot
%------------------------------------------------------------------------
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

% Preallocate outputs
u=nan(size(igood,1),nLayers,length(plotOPTS.Time_record));
v=nan(size(igood,1),nLayers,length(plotOPTS.Time_record));

% Check if we're running
for aa=1:length(plotOPTS.Time_record)
    fprintf('Time step %i of %i\n',aa,length(plotOPTS.Time_record))


    % Mesh goes underneath vectors.
    if plotOPTS.do_mesh
        % plot vertices
        [X,Y]=m_ll2xy(plotOPTS.mesh.lon,plotOPTS.mesh.lat,'clip','on');
        Plots(plotOPTS.figure).mesh=patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
            'EdgeColor',[0.6 0.6 0.6],'FaceColor','none');hold on
    end

    for ii=1:nLayers
        u(:,ii)=squeeze(FVCOM.u(igood,nLayersRange(ii),(aa)));
        v(:,ii)=squeeze(FVCOM.v(igood,nLayersRange(ii),(aa)));
    end
    for jj=1:nLayers
        [Plots(plotOPTS.figure).handles(jj),~]=m_vec(plotOPTS.vel_sca,x(igood),y(igood),...
            u(:,jj),v(:,jj),colourSpec(jj,:),...
            'shaftwidth',1,'headwidth',2);
        hold on
    end
    pause(plotOPTS.pause)
    if aa~=length(plotOPTS.Time_record)
        delete(Plots(plotOPTS.figure).handles(:))
    end
end
