function [Plots]=do_ph_vertical_profile(plotOPTS,FVCOM,transectPoints)
% For a supplied transect, plot the vertical profile through the water
% column supplied in FVCOM.var_plot.

% Get the relevant time intervals
if isempty(plotOPTS.save_intervals)
    plotOPTS.save_intervals=1:length(plotOPTS.Time_record);
    dontSave=1;
else
    dontSave=0;
end

if isfield(plotOPTS,'altColours') && plotOPTS.altColours==1
    colourSpec=flipud(jet);
else
    clear nColours nColourIn nColourOut colourSpec
    % Build a colour palette which matches Jerry's ranges.
    % Dark Red -> Red -> Amber -> Yellow -> Green -> Blue
    nColours=200;
    nColourIn=[1,nColours*0.15,nColours*0.6,nColours*0.75,nColours*0.9,nColours];
    nColourOut=1:nColours; % Gives a nice continuous palette.
    %                          DR    R     A    Y    G     B
    cRed=interp1(nColourIn,  [0.62, 0.9, 1    , 1  , 0  , 0.46],nColourOut);
    cGreen=interp1(nColourIn,[0   , 0  , 0.52 , 1  , 0.8, 0.63],nColourOut);
    cBlue=interp1(nColourIn, [0.2 , 0.2, 0    , 0  , 0  , 0.83],nColourOut);
    colourSpec=[cRed;cGreen;cBlue]';
end

dataToUse=single(FVCOM.(plotOPTS.var_plot));

saveInc=1;
for aa=1:length(plotOPTS.Time_record)

    verticalProfile=squeeze(dataToUse(transectPoints.trn_nodes.idx,:,plotOPTS.save_intervals(aa)));
    
%     xValues=repmat(transectPoints.trn_dis,1,size(verticalProfile,2));
    % Not sure about those distance values...
    xNorm=transectPoints.trn_nodes.x-min(transectPoints.trn_nodes.x);
    yNorm=transectPoints.trn_nodes.y-min(transectPoints.trn_nodes.y);
    xDist=sqrt(xNorm.^2+yNorm.^2);
    xValues=repmat(xDist,1,size(verticalProfile,2));

    zeta=squeeze(FVCOM.zeta(transectPoints.trn_nodes.idx,plotOPTS.save_intervals(aa)));
    waterDepth=FVCOM.h(transectPoints.trn_nodes.idx);
    zValues=(zeta+waterDepth)*FVCOM.siglay(1,:);

    % plot profile
    fprintf('Time step %i of %i\n',aa,length(plotOPTS.Time_record))
    figure(plotOPTS.figure); clf
    hold on
    Plots(1).handles=contourf(...
        xValues/1000,... % km
        zValues,...
        fliplr(verticalProfile),... % get vertical right way up
        200,'edgecolor','none');
%     ylim([-max(FVCOM.h(:)) -0.5])
    xlabel('Distance along transect (km)')
    ylabel('Depth (m)')
    colormap(colourSpec)
    ch=colorbar;
    set(ch,'FontSize',10);
    ylabel(ch,'pH')

    if ~isempty(plotOPTS.vlims)
        caxis(plotOPTS.vlims);
        % Get sensible tick formatting (if necessary)
        if plotOPTS.vlims(2)-plotOPTS.vlims(1)<1e-4
            xx=0:length(plotOPTS.vlims)-1;
            yy=plotOPTS.vlims;
            xxi=0:1/6:1; % six ticks
            cticks=interp1(xx,yy,xxi);
            set(ch,'YTick',double(cticks))
            set(ch,'yticklabel',cticks)
        end
    end



    pause(plotOPTS.pause)
    if plotOPTS.save_output % Are we even trying to save figures?
        if saveInc<=length(plotOPTS.save_intervals) && plotOPTS.save_intervals(aa)==plotOPTS.save_intervals(saveInc) && ~dontSave
            % Save output
            fprintf('Saving figure... ')
            set(findobj(gcf,'Type','text'),'FontSize',10)
            %set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize'))) 
            set(gcf,'PaperPositionMode','auto');
            set(gcf,'renderer','painters'); % for vector output in pdfs
            print(gcf,'-dpdf','-r600',[plotOPTS.FVCOM_plot_dir,plotOPTS.var_plot,'/pdf/',plotOPTS.fig_name,'_vertical_profile_',plotOPTS.var_plot,'_',num2str(plotOPTS.save_intervals(aa)),'.pdf']); % pdf
            %print(gcf,'-dpng','-r600',[plotOPTS.FVCOM_plot_dir,plotOPTS.var_plot,'/png/',plotOPTS.fig_name,'_layer=',num2str(plotOPTS.nz_plot),'_',plotOPTS.var_plot,'_change_',num2str(plotOPTS.save_intervals(aa)),'.png']); % png
            saveInc=saveInc+1;
            fprintf('done.\n')
        end
    end

    if aa~=length(plotOPTS.Time_record)
        clf
    end
end

return