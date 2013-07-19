function [Plots]=do_ph_max_change_plot(plotOPTS,FVCOM,startIdx)
% Calculate the change in pH and plot accordingly.
m_mappath;

% Check we have some of the required fields.
if ~isfield(FVCOM,plotOPTS.var_plot)
    error('Need %s input to calculate change in %s.',plotOPTS.var_plot,plotOPTS.var_plot)
end

% Build a colour palette which matches Jerry's ranges.
% Dark Red -> Red -> Amber -> Yellow -> Green -> Blue
nColours=200;
nColourIn=[1,nColours*0.15,nColours*0.6,nColours*0.75,nColours*0.9,nColours];
nColourOut=1:nColours; % Gives a nice continuous palette.
%                          DR    R     A    Y    G     B
cRed=interp1(nColourIn,  [0.62, 0.9, 1    , 1  , 0  , 0.46],nColourOut);
cGreen=interp1(nColourIn,[0   , 0  , 0.52 , 1  , 0.8, 0.63],nColourOut);
cBlue=interp1(nColourIn, [0.2 , 0.2, 0    , 0  , 0  , 0.83],nColourOut);
colourSpec=flipud([cRed;cGreen;cBlue]'); % flip since we have negative scale

if isfield(plotOPTS,'altColours') && plotOPTS.altColours==1
    clear nColours nColourIn nColourOut colourSpec
    % Build a colour palette which matches Jerry's ranges.
    % Dark Red -> Red -> Amber -> Yellow -> Green -> Blue
    nColours=200;
    nColourIn=[1,nColours*0.1,nColours*0.2,nColours*0.3,nColours*0.4,nColours];
    nColourOut=1:nColours; % Gives a nice continuous palette.
    %                          DB    B     LB   DG    G     LG
    cRed=interp1(nColourIn,  [0   , 0   , 0   , 0.04, 0 , 0.2 ],nColourOut);
    cGreen=interp1(nColourIn,[0   , 0   , 0.2 , 0.51, 0.6, 0.76 ],nColourOut);
    cBlue=interp1(nColourIn, [0.4 , 0.6 , 0.79, 0.78, 0.6  , 0],nColourOut);
    colourSpec=[cRed;cGreen;cBlue]';
end

dataToUse=single(FVCOM.(plotOPTS.var_plot));

% Get the background condition. Depth average if necessary.
if isfield(plotOPTS,'depth_average') && plotOPTS.depth_average
    bgPH=squeeze(mean(dataToUse(:,:,startIdx),2));
else
    bgPH=squeeze(dataToUse(:,plotOPTS.nz_plot,startIdx));
end

% Are we depth averaging?
if isfield(plotOPTS,'depth_average') && plotOPTS.depth_average
    phData=squeeze(mean(dataToUse,2));
else
    phData=squeeze(dataToUse(:,plotOPTS.nz_plot,:));
end

phDiff=phData-repmat(bgPH,1,size(phData,2));

% Are we depth averaging?
% if isfield(plotOPTS,'depth_average') && plotOPTS.depth_average
%     % Calculate the difference for each time step from the initial
%     % conditions.
%     phDiff=squeeze(mean(FVCOM.(plotOPTS.var_plot),2))-repmat(bgPH,1,length(plotOPTS.Time_record));
% else
%     phDiff=squeeze(FVCOM.(plotOPTS.var_plot)(:,plotOPTS.nz_plot,:))-repmat(bgPH,1,size(FVCOM.(plotOPTS.var_plot),3));
% end

% Not doing cumulative differences any more
% Check if we're doing the cumulative difference
% if isfield(plotOPTS,'summed_ph') && plotOPTS.summed_ph
%     phDiff=cumsum(phDiff,2);
% end

% Do the maximum change from initial conditions across all time steps (i.e.
% not the cumulative change). abs() the difference to find the greatest
% change, either increase or decrease.
% We're also averaging here over the number of timesteps so that we can fix
% the colour palette to something sensible.
phDiffMax=min((phDiff),[],2); %./(size(plotOPTS.Time_record,2)*(plotOPTS.Time_record(2)-plotOPTS.Time_record(1)));
figure(plotOPTS.figure); clf
m_proj('UTM','lon',[plotOPTS.range_lon],'lat',[plotOPTS.range_lat],'zon',plotOPTS.zone,'ell','grs80')
m_grid('box','fancy')
m_usercoast(plotOPTS.coastline_file,'Color','k','LineWidth',3);
[X,Y]=m_ll2xy(plotOPTS.mesh.lon,plotOPTS.mesh.lat,'clip','on');

% plot map with pH change
hold on
Plots(1).handles=patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
    'Cdata',phDiffMax,'edgecolor','interp','facecolor','interp');
caxis(plotOPTS.clims)
colormap(flipud(colourSpec))
ch=colorbar;
set(ch,'FontSize',10);
ylabel(ch,'pH change');
% ylabel(ch,[plotOPTS.var_plot,' change'])
% check if mesh elements are required
if plotOPTS.do_mesh
    % plot vertices
    [X,Y]=m_ll2xy(plotOPTS.mesh.lon,plotOPTS.mesh.lat,'clip','on');
    patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
        'EdgeColor',[0.6 0.6 0.6],'FaceColor','none'); hold on
end

if plotOPTS.save_output % Are we even trying to save figures?
    % Save output
    fprintf('Saving figure... ')
    set(findobj(gcf,'Type','text'),'FontSize',10)
    %set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')))
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'renderer','painters'); % for vector output in pdfs
    print(gcf,'-dpdf','-r600',[plotOPTS.FVCOM_plot_dir,plotOPTS.var_plot,'/pdf/',plotOPTS.fig_name,'_layer=',num2str(plotOPTS.nz_plot),'_',plotOPTS.var_plot,'_max_change.pdf']); % pdf
    %print(gcf,'-dpng','-r600',[plotOPTS.FVCOM_plot_dir,plotOPTS.var_plot,'/png/',plotOPTS.fig_name,'_layer=',num2str(plotOPTS.nz_plot),'_',plotOPTS.var_plot,'_max_change.png']); % png
    fprintf('done.\n')
end


return
