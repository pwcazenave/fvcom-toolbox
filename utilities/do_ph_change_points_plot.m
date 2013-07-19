function [Plots]=do_ph_change_points_plot(plotOPTS,FVCOM,startIdx)
% Calculate the change in the given parameter (plotOPTS.var_plot) and plot
% accordingly.
m_mappath;

plotOPTS.figure=1;
plotOPTS.var_plot='ph';
plotOPTS.do_mesh=1; % Add the mesh?
startIdx=1; %:size(FVCOM.zeta,2);

plotOPTS.threshold_change=-0.2;
% near bottom, but not actually the bottom since it doesn't change at all
plotOPTS.nz_plot=1; 

[nx,nz,ttot]=size(FVCOM.(plotOPTS.var_plot));

% Check we have some of the required fields.
if ~isfield(FVCOM,plotOPTS.var_plot)
    error('Need %s input to calculate change in %s.',plotOPTS.var_plot,plotOPTS.var_plot)
end

% Get the relevant time intervals
if isfield(plotOPTS,'save_intervals') && isempty(plotOPTS.save_intervals)
    plotOPTS.save_intervals=1:length(plotOPTS.Time_record);
    dontSave=1;
else
    dontSave=0;
end

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
        % I'm pretty certain this needs to integrate for the time between
        % time steps, otherwise halving the output time step would halve
        % the exposed volume, which can't be right...
        % However, if we know that the threshold has been triggered n
        % times, we just need n times the volume, irrespective of time.
        % Right? 
        cellVolume(:,ii,jj)=totalDepth(:,jj).*sigThickness(:,ii).*FVCOM.art1; %*dt;
    end
end

% Calculate the difference between the current step and the previous step.
% phChange=diff(squeeze(FVCOM.(plotOPTS.var_plot)(:,plotOPTS.nz_plot,:)),[],2);
% Nope, get the background condition and compare against that. Depth
% average if necessary.
if isfield(plotOPTS,'depth_average') && plotOPTS.depth_average
    bgPH=squeeze(mean(FVCOM.(plotOPTS.var_plot)(:,:,startIdx),2));
else
    bgPH=squeeze(FVCOM.(plotOPTS.var_plot)(:,plotOPTS.nz_plot,startIdx));
end

if isfield(plotOPTS,'depth_average') && plotOPTS.depth_average
    % Now, for each successive time step, calculate the difference between the
    % current time step and the background level.
    phMean=squeeze(mean(FVCOM.(plotOPTS.var_plot),2));
    phDiff=phMean-repmat(bgPH,1,size(phMean,2));
else
    phDiff=squeeze(FVCOM.(plotOPTS.var_plot)(:,plotOPTS.nz_plot,:))-repmat(bgPH,1,size(FVCOM.(plotOPTS.var_plot),3));
end

% Try a cumulative difference
phDiffCumulative=zeros(nx,ttot,'single');
for tt=1:ttot
    if tt>1 && tt<ttot
        currpH=FVCOM.(plotOPTS.var_plot)(:,plotOPTS.nz_plot,tt);
        nextpH=FVCOM.(plotOPTS.var_plot)(:,plotOPTS.nz_plot,tt+1);
        phDiffCumulative(:,tt)=phDiffCumulative(:,tt)+(currpH-nextpH);
    elseif tt==1
        currpH=FVCOM.(plotOPTS.var_plot)(:,plotOPTS.nz_plot,tt);
        nextpH=FVCOM.(plotOPTS.var_plot)(:,plotOPTS.nz_plot,tt+1);
        phDiffCumulative(:,tt)=currpH-nextpH;
    end
end

[X,Y]=m_ll2xy(plotOPTS.mesh.lon,plotOPTS.mesh.lat,'clip','on');

foundXY=cell(1,length(plotOPTS.nz_plot));
countSites=zeros(nx,1);

for tt=1:ttot
    % Find drops greater than the threshold
    idx=find(phDiff(:,tt)<plotOPTS.threshold_change);
    if ~isempty(idx)
        countSites(idx)=countSites(idx)+1;
        foundXY{tt}.xy=[X(idx),Y(idx)];
    end
end

% Now we can use countSites to calculate the volume exposed to the change,
% although this ignores the effect of the tide (because we don't know which
% time step triggered the threshold condition test)
totalVolume=squeeze(cellVolume(:,plotOPTS.nz_plot,1)).*countSites;

% Visualise the count threshold
Plots(1).handles=patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
    'Cdata',countSites,...
    'edgecolor','interp','facecolor','interp');
if plotOPTS.do_mesh
    % plot vertices
    [X,Y]=m_ll2xy(plotOPTS.mesh.lon,plotOPTS.mesh.lat,'clip','on');
    patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
        'EdgeColor',[0.6 0.6 0.6],'FaceColor','none'); hold on
end
colorbar
    
figure;
patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
    'Cdata',totalVolume/1e9,...
    'edgecolor','interp','facecolor','interp');
if plotOPTS.do_mesh
    % plot vertices
    [X,Y]=m_ll2xy(plotOPTS.mesh.lon,plotOPTS.mesh.lat,'clip','on');
    patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
        'EdgeColor',[0.6 0.6 0.6],'FaceColor','none'); hold on
end
colorbar
    
figure;
patch('Vertices',[X,Y],'Faces',plotOPTS.mesh.tri,...
    'Cdata',squeeze(FVCOM.ph(:,1,2)),...
    'edgecolor','interp','facecolor','interp');
colorbar

% figure(1)
% subplot(2,1,1)
% plot(Time_record-min(Time_record),squeeze(FVCOM.ph(1316,plotOPTS.nz_plot,:)))
% subplot(2,1,2)
% plot(Time_record(1:end-1)-min(Time_record),phChange(1316,:),'r')
% hold on
% plot(Time_record(idxAtLeak)-min(Time_record),phChange(1316,idxAtLeak),'g.')
% % Add the threshold
% plot([min(Time_record)-min(Time_record),max(Time_record)-min(Time_record)],[plotOPTS.threshold_change,plotOPTS.threshold_change],'k--')
