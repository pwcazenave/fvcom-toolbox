function [Plots]=do_vector_plot(plotOPTS,FVCOM)
% Plot vectors of parameters of interest.
m_mappath;

figure(plotOPTS.figure); clf
m_proj('UTM','lon',[plotOPTS.range_lon],'lat',[plotOPTS.range_lat],'zon',30,'ell','grs80')
m_grid('box','fancy')
m_usercoast(plotOPTS.coastline_file,'Color','k','LineWidth',3);
[x,y]=m_ll2ll(FVCOM.xc,FVCOM.yc); x=x+6;

igood = find (x < plotOPTS.range_lon(2) & x > plotOPTS.range_lon(1) &...
    y < plotOPTS.range_lat(2) & y > plotOPTS.range_lat(1));
igood=igood(1:plotOPTS.data_dec:end);

if isfield(plotOPTS,'nz_plot_vec')
    nLayers=size(plotOPTS.nz_plot_vec,2);
    nLayersRange=plotOPTS.nz_plot_vec;
else
    nLayers=size(plotOPTS.nz_plot,2);
    nLayersRange=plotOPTS.nz_plot;
end

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
        Plots(plotOPTS.figure).handles=patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
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
