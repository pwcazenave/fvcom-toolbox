function [Plots,totalVol]=do_volume(plotOPTS,FVCOM,startIdx,thresholdValue)
% Calculate a volume for a given time step above a threshold value.
m_mappath;

[nx,nz,ttot]=size(FVCOM.(plotOPTS.var_plot));

% Get the time step (in seconds)
dt=round((plotOPTS.Time_record(2)-plotOPTS.Time_record(1))*24*60*60);

% Sigma layer fraction (nz*nz)
sigThickness=roundn(abs(diff(FVCOM.siglev,1,2)),-5); % roundn to even values out.
% Total depth (nz*ttot)
totalDepth=repmat(FVCOM.h,1,size(FVCOM.zeta,2))+FVCOM.zeta;
% Volume for every element for each time step
cellVolume=nan(nx,nz,ttot);
for ii=1:nz
    for jj=1:ttot
        cellVolume(:,ii,jj)=totalDepth(:,jj).*sigThickness(:,ii).*FVCOM.art1*dt;
    end
end

if plotOPTS.nz_plot==0
    % Depth averaging the results
    cellVolume=squeeze(sum(cellVolume,2));
    dataArray=squeeze(mean(FVCOM.(plotOPTS.var_plot)(:,:,:),2));
else
    dataArray=FVCOM.(plotOPTS.var_plot)(:,:,:);
end

%%
% Calculate volume of elements matching threshold condition
totalVolTemp=0;
colourSpec=hsv(length(plotOPTS.nz_plot));
plotSymbols={'.','o','x','+','^','*','p','h'};

for tt=1:size(startIdx,2)
    
    foundXY=cell(1,length(plotOPTS.nz_plot));
    
    figure(plotOPTS.figure); clf
    m_proj('UTM','lon',[plotOPTS.range_lon],'lat',[plotOPTS.range_lat],'zon',plotOPTS.zone,'ell','grs80')
    m_grid('box','fancy')
    m_usercoast(plotOPTS.coastline_file,'Color','k','LineWidth',3);
    [X,Y]=m_ll2xy(plotOPTS.mesh.lon,plotOPTS.mesh.lat,'clip','on');

    fprintf('Time step %i of %i\n',startIdx(tt),length(plotOPTS.Time_record))
    % Find all values of pH greater than the thresholdValue for the layers
    % specified in .nz_plot only.
    
    for ss=1:size(plotOPTS.nz_plot,2)
        if plotOPTS.nz_plot==0
            % Depth averaged, so easy.
            idx=thresholdTest(dataArray(:,startIdx(tt)),plotOPTS.volume_thresh,thresholdValue);
            % Do we have any locations where the condition has been met?
            if ~isempty(idx)
                totalVolTemp=totalVolTemp+sum(sum(cellVolume(idx,startIdx(tt))));
            end
        else
            % Loop through the layers
            for qq=1:size(plotOPTS.nz_plot,2)
                idx=thresholdTest(dataArray(:,plotOPTS.nz_plot(qq),startIdx(tt)),plotOPTS.volume_thresh,thresholdValue);
                % Again, look for areas where we have met the threshold
                % condition. Here, we need to save the xy coordinates for
                % each find with the layer number.
                if ~isempty(idx)
                    totalVolTemp=totalVolTemp+sum(sum(cellVolume(idx,plotOPTS.nz_plot(qq),startIdx(tt))));
                    % This may skip the odd layer...
                    foundXY{qq}.xy=[X(idx),Y(idx),repmat(plotOPTS.nz_plot(qq),length(idx),1)];
                end
            end
        end
    end
    if tt==1;
        totalVol(tt)=totalVolTemp;
    else
        totalVol(tt)=totalVol(tt-1)+totalVolTemp;
    end
    % Do a figure identifying the locations above the threshold for the
    % surface sigma layer.
    if plotOPTS.nz_plot==0
        Plots(plotOPTS.figure).handles=patch('Vertices',[X,Y],'Faces',...
            plotOPTS.mesh.tri,'Cdata',squeeze(cellVolume(:,tt)),...
            'edgecolor','interp','facecolor','interp');
    else
            Plots(plotOPTS.figure).handles=patch('Vertices',[X,Y],'Faces',...
        plotOPTS.mesh.tri,'Cdata',squeeze(cellVolume(:,1,tt)),...
        'edgecolor','interp','facecolor','interp');
    end
    hold on
    % TODO: For multiple layers, do unique indices for each layer so they
    % can be plotted separately here (a la do_vector_plot.m).
    if ~isempty(idx)
        for kk=1:size(plotOPTS.nz_plot,2);
            plot(foundXY{kk}.xy(:,1),foundXY{kk}.xy(:,2),plotSymbols{mod(kk,length(plotSymbols)+1)},'Color',colourSpec(kk,:))
        end
    end
    if plotOPTS.do_mesh
        % plot vertices
        [X,Y]=m_ll2xy(plotOPTS.mesh.lon,plotOPTS.mesh.lat,'clip','on');
        patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
            'EdgeColor',[0.6 0.6 0.6],'FaceColor','none'); hold on
    end
    % Some useful text
    [textX,textY]=m_ll2xy(-4.45,50.13);
    text(textX,textY,sprintf('Total volume:\t\t%.2fkm^{3}\nVolume this time step:\t%.2fkm^{3}\n',totalVol(tt)/1e9,totalVolTemp/1e9))
        
    pause(plotOPTS.pause)
    caxis(plotOPTS.clims)
    colorbar
    set(get(colorbar,'YLabel'),'String','Volume (m^{3})')
    
    if tt~=size(startIdx,2)
        delete(Plots(plotOPTS.figure).handles)
    end
end

function idx=thresholdTest(dataArray,thresholdType,thresholdValue)
% Abstract some of the complexity into a function. 
% Usage: thresholdTest(inputData,thresholdValue)
    switch thresholdType
        case -1
            idx=find(dataArray<thresholdValue);
        case 1
            idx=find(dataArray>thresholdValue);
        case 0
            idx=find(dataArray==thresholdValue);
        otherwise
            error('Unrecognised value for ''plotOPTS.volume_thresh''.')
    end

    return