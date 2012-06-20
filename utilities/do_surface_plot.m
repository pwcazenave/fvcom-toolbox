function [Plots]=do_surface_plot(plotOPTS,FVCOM)
% reads image and plots tracks or stations
m_mappath;

figure(plotOPTS.figure); clf
m_proj('UTM','lon',[plotOPTS.range_lon],'lat',[plotOPTS.range_lat],'zon',plotOPTS.zone,'ell','grs80')
m_grid('box','fancy')
%m_usercoast(plotOPTS.coastline_file,'Color','k','LineWidth',3);
%[X,Y]=m_ll2xy(plotOPTS.mesh.lon,plotOPTS.mesh.lat,'clip','on');
X=plotOPTS.mesh.x;
Y=plotOPTS.mesh.y;
for aa=1:length(plotOPTS.Time_record)
    % plot map with plotOPTS.var_plot
    hold on
    try
        % 3D data (i.e. DYE, pH etc.)
        Plots(1).handles=patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
            'Cdata',squeeze(FVCOM.(plotOPTS.var_plot)(:,plotOPTS.nz_plot,aa)),...
            'edgecolor','interp','facecolor','interp');
    catch
        % 2D data only (i.e. zeta etc.)
        Plots(1).handles=patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
            'Cdata',squeeze(FVCOM.(plotOPTS.var_plot)(:,aa)),...
            'edgecolor','interp','facecolor','interp');
    end
    fprintf('Time step %i of %i\n',aa,length(plotOPTS.Time_record))
    caxis(plotOPTS.clims)
    colorbar
    set(get(colorbar,'YLabel'),'String',plotOPTS.var_plot)
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
    pause(plotOPTS.pause)
    if aa~=length(plotOPTS.Time_record)
        delete(Plots(1).handles)
    end
end

return

