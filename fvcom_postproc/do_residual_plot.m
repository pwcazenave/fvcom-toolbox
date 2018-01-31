function [Plots]=do_residual_plot(plotOPTS,FVCOM,dt)
% Take the output of do_residual and plot as a vector figure. Summarises a
% specified interval of time as a single long-term direction and magnitude.
m_mappath;

warning('on','FVCOM:Plot:ResidualAnalysis')

figure(plotOPTS.figure); clf
m_proj('UTM','lon',[plotOPTS.range_lon],'lat',[plotOPTS.range_lat],'zon',30,'ell','grs80');
m_grid('box','fancy');
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

% Check we're not depth averaging values.
if isfield(plotOPTS,'depth_average') && plotOPTS.depth_average
    if nLayers>1
        warning('FVCOM:Plot:ResidualAnalysis','Depth averaging has been set as well as a specific number of layers to extract. Usually one or the other is preferred.')
    end
    nLayers=1;
    nLayersRange=1;
    % Average through all depths. Don't squeeze() here as that's taken
    % care of in do_residual().
    uIn=mean(FVCOM.u,2);
    vIn=mean(FVCOM.v,2);
else
    uIn=FVCOM.u(:,nLayersRange,:);
    vIn=FVCOM.v(:,nLayersRange,:);
end

if nLayers==1
    colourSpec=[0 0 0];
else
    colourSpec=colormap(hsv(nLayers));
    setColourMap=1;
end

% We're not using uRes and vRes here, but if you wanted to do a PVD, then
% you would use:
%   plot(uRes(someElement,someLayer,:),vRes(someElement,someLayer,:),'.-'),
% for example.
[rDir,rMag,uRes,vRes]=do_residual(uIn,vIn,dt);

% Mesh goes underneath the vectors
if plotOPTS.do_mesh
    % plot vertices
    [X,Y]=m_ll2xy(plotOPTS.mesh.lon,plotOPTS.mesh.lat,'clip','on');
    Plots(plotOPTS.figure).handles=patch('Vertices',[X,Y],...
        'Faces',plotOPTS.mesh.tri,'EdgeColor',[0.6 0.6 0.6],...
        'FaceColor','none'); hold on
end

for ii=1:nLayers
    % Decompose to vector components for m_vec.
    uVec=rMag(:,ii,:).*sind(rDir(:,ii,:));
    vVec=rMag(:,ii,:).*cosd(rDir(:,ii,:));
    [Plots(plotOPTS.figure).handles(ii),~]=m_vec(plotOPTS.vel_sca,...
        x(igood),y(igood),squeeze(uVec(igood)),squeeze(vVec(igood)),...
        colourSpec(ii,:),'shaftwidth',1,'headwidth',2);
    if exist('setColourMap','var')
        if setColourMap
            colorbar
            set(get(colorbar,'YLabel'),'String','Layer')
        end
    end
end

