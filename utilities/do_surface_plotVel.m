function [Plots]=do_surface_plotVel(plotOPTS,FVCOM)
% reads image and plots tracks or stations
m_mappath;

figure(plotOPTS.figure)
% set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5])
clf
m_proj('UTM','lon',[plotOPTS.range_lon],'lat',[plotOPTS.range_lat],'zon',plotOPTS.zone,'ell','grs80')
m_grid('box','fancy')
set(gca,'FontSize',18)

m_usercoast(plotOPTS.coastline_file,'Color','k','LineWidth',3);
[X,Y]=m_ll2xy(plotOPTS.mesh.lon,plotOPTS.mesh.lat,'clip','on');


if isfield(plotOPTS,'Time_record') & isfield(FVCOM,'mattime')
    if length(plotOPTS.Time_record)==1
        [dump,igoodT] = min(abs(plotOPTS.Time_record(1)-FVCOM.mattime));
    else
        
        igoodT = find( plotOPTS.Time_record(1) <= FVCOM.mattime &   FVCOM.mattime <=  plotOPTS.Time_record(end) );
    end
else
    igoodT=(1:length(FVCOM.(plotOPTS.var_plot)));
end


for aa=igoodT
    % plot map with plotOPTS.var_plot
    hold on
    try
        % 3D data (i.e. DYE, pH etc.)
        Plots(1).handles=patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
            'FaceVertexCdata',squeeze(FVCOM.(plotOPTS.var_plot)(:,plotOPTS.nz_plot,aa)),...
            'CDataMapping','scaled','edgecolor','none','FaceColor','flat');
%          Plots(1).handles=scatter(plotOPTS.mesh.xc,plotOPTS.mesh.yc,20,squeeze(FVCOM.(plotOPTS.var_plot)(:,plotOPTS.nz_plot,aa)),'filled')
    catch
        % 2D data only (i.e. zeta etc.)
        Plots(1).handles=patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
            'FaceVertexCdata',squeeze(FVCOM.(plotOPTS.var_plot)(:,aa)),...
            'CDataMapping','scaled','edgecolor','none','FaceColor','flat');
%          Plots(1).handles=scatter(plotOPTS.mesh.xc,plotOPTS.mesh.yc,20,squeeze(FVCOM.(plotOPTS.var_plot)(:,aa)),'filled')
    end
    fprintf('Time step %i of %i\n',aa,length(plotOPTS.Time_record))
    display(['Time ',datestr(FVCOM.mattime(aa))])
    caxis(plotOPTS.clims)
    colorbar
    set(get(colorbar,'YLabel'),'String',plotOPTS.var_plot,'FontSize',14)
    % check if mesh elements are required
    if plotOPTS.do_mesh
        % plot vertices
        [X,Y]=m_ll2xy(plotOPTS.mesh.lon,plotOPTS.mesh.lat,'clip','on');
        patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
            'EdgeColor',[0.6 0.6 0.6],'FaceColor','none'); hold on
    end
    %----------------------------------------------------------------------
    % Only in my case it needs adding 6 because proj automatically
    % determines a reference long in strides of 6deg while m_map doesn't
    %----------------------------------------------------------------------
   m_usercoast(plotOPTS.coastline_file,'patch',[0.6, 0.6, 0.6]);
   pause(plotOPTS.pause)
%     if aa~=length(plotOPTS.Time_record)
%         delete(Plots(1).handles)
%     end
end

return

