function [Plots]=do_transect_plot(plotOPTS,FVCOM)
%
% Function to display transect plots of FVCOM variables (i.e. temperature)
%
%  [Plots]=do_transect_plot(plotOPTS,FVCOM)
%
% DESCRIPTION:
%    Generates transect plots of variables
%
% INPUT:
%   plotOPTS   = structure array with predefined options for generating the
%   maps
%   FVCOM  = data extracted from FVCOM NC file. See read_netCDF_FVCOM for
%   details
%
%   plotOPTS.range_lat: [50.1000 50.4000]
%   plotOPTS.range_lon: [-4.5000 -3.8500]
%   If the above exist, it will generate an m_map plot with transect
%   position
%   plotOPTS.mesh: [1x1 struct]
%   plotOPTS.coastline_file: '../mat/tamar3_0coast.mat'
%   plotOPTS.zone: 30
%   plotOPTS.ell: 'grs80'
%   plotOPTS.fig_name: 'co2_S5_slowleak'
%   plotOPTS.mesh: [1x1 struct]
%   plotOPTS.var_plot: 'PH'
%   plotOPTS.clims: [6 8]
%   plotOPTS.nz_plot: 1
%   plotOPTS.figure: 1
%   plotOPTS.Time_record: 7.3271e+05
%   plotOPTS.data_dec: 5
%   plotOPTS.vel_sca: 5
%   plotOPTS.pause: 0.5000

%
% OUTPUT:
%   Plots = structure array with figure handles
%
% EXAMPLE USAGE
%    [Plots]=do_transect_plot(plotOPTS,FVCOM)
%
% Author(s):
%    Ricardo Torres (Plymouth Marine Laboratory)
%
% Revision history
%
%==============================================================================
figure(plotOPTS.figure);
if ~(isfield(FVCOM,'zeta') | ~isfield(FVCOM,'h') | ~isfield(FVCOM,'siglay'))
    fprintf('You seem to be missing either  zeta, h, or siglay from FVCOM variable\n');
    fprintf(' Stoping. Try again\n');
    Plots=[];
    return
end

% convert sigma levels to z levels
nz=size(FVCOM.(plotOPTS.var_plot),2)
if isfield(plotOPTS,'Time_record')
    nt=length(plotOPTS.Time_record)
else
    nt=1;
end

nxy=length(plotOPTS.trn_nodes.y)
FVCOM.z=ones(nxy,nz,nt);
for tt=1:nt
    z=FVCOM.h(:)+FVCOM.zeta(:,tt);
    for zz=1:nz
        FVCOM.z(:,zz,tt)=z+z.*FVCOM.siglay(:,zz);
    end
end
% plot animation of timeseries of transect
XX=repmat(plotOPTS.trn_dis(:)',nz,1)';
for tt=1:nt
    YY=squeeze(FVCOM.z(:,:,tt));
    ZZ=squeeze(FVCOM.(plotOPTS.var_plot)(:,:,tt));
    Plots(plotOPTS.figure).handles=contourf(XX/100000,-YY,(ZZ),200,'edgecolor','none');
    caxis(plotOPTS.clims);
    colorbar
    xlabel('Distance along Transect')
    ylabel('Depth m')
    ylim([-max(FVCOM.h(:)) 0])
    pause(.5)
end


if isfield(plotOPTS,'range_lat')
    m_mappath;
    %plot map with transect overlaid
    figure(plotOPTS.figure+1);clf
    m_proj('UTM','lon',[ plotOPTS.range_lon],'lat',[plotOPTS.range_lat],'zon',plotOPTS.zone,'ell',plotOPTS.ell)
    m_grid('box','fancy')
    % add coastline if present
    if (isfield(plotOPTS,'coastline_file') && ~isempty(plotOPTS.coastline_file) )
        m_usercoast(plotOPTS.coastline_file,'Color','k','LineWidth',3);
    end
    
    %plot vertices
    col=0.7*[1 1 1];
    if plotOPTS.do_mesh
        %plot vertices
        [X,Y]=m_ll2xy(plotOPTS.mesh.lon,plotOPTS.mesh.lat,'clip','on');
        Plots(plotOPTS.figure).mesh=patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
            'EdgeColor',[0.6 0.6 0.6],'FaceColor','none');hold on
        
    end
    % plot with matlab
    hold on
    [x,y]=m_ll2ll(plotOPTS.trn_nodes.x,plotOPTS.trn_nodes.y);
    [X,Y]=m_ll2xy(x+6,y,'clip','on');
    plot(X,Y,'r','LineWidth',2.5)
end

