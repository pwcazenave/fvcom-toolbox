function [totalVolume] = do_volume_change(plotOPTS,FVCOM)
% DO_VOLUME_CHANGE Calculate volume of water which experiences a change in
%   pH beyond some defined threshold.

% Make a 3D array of the volumes to get to total volume which experiences
% the change in pH.
[nx,nz,ttot]=size(FVCOM.(plotOPTS.var_plot));

% Get the time step (in seconds)
dt=round((plotOPTS.Time_record(2)-plotOPTS.Time_record(1))*24*60*60);
% Make sure we're actually able to sample at the rate requested.
if plotOPTS.change_type~=0 && dt/(60*60)>plotOPTS.change_type
    error('Output file sampling frequency is coarser than specified time sampling.')
end

% Sigma layer fraction (nz*nz)
sigThickness=roundn(abs(diff(FVCOM.siglev,1,2)),-5); % roundn to even values out.
% Total depth (nz*ttot)
totalDepth=repmat(FVCOM.h,1,size(FVCOM.zeta,2))+FVCOM.zeta;
% Volume for every element for each time step
cellVolume=nan(nx,nz,ttot);
for ii=1:nz
    for jj=1:ttot
        cellVolume(:,ii,jj)=totalDepth(:,jj).*sigThickness(:,ii).*FVCOM.art1.*dt;
    end
end

if plotOPTS.change_type==0 || plotOPTS.change_type==dt/(60*60)
    % Instantaneous change (i.e. a change between outputs beyond a given value.
    phChange=diff(FVCOM.DYE,1,3);
    
    % Find the elements with the change beyond the threshold.
    idx=phChange<plotOPTS.threshold_change;

    % Using those indices, find the volume which is affected.
    totalVolume=sum(cellVolume(idx));
else
    % The more complicated n-hour change.
    
    % Use the time period over which we're interested in to get an index
    % jump value.
    dtJump=(plotOPTS.change_type*60*60)/dt;
    if dtJump-round(dtJump)~=0
        error('Output sampling is not compatible with time period of change.')
    end

    phChange=nan(nx,nz,ttot-dtJump);
    for tt=1:length(plotOPTS.Time_record)
        if tt<=length(plotOPTS.Time_record)-dtJump
            phChange(:,:,tt)=(FVCOM.DYE(:,:,tt+dtJump)-FVCOM.DYE(:,:,tt));
        end
    end
    
    % Find the change values that happen for at least n-hours.
    % Probably best to look away now if you're interested in clean quick
    % code.
    totalVolume=0;
    for ii=1:nx
        for jj=1:nz
            totalVolumeCurrent=get_runs(plotOPTS,squeeze(phChange(ii,jj,:)),squeeze(cellVolume(ii,jj,:)),plotOPTS.change_type,plotOPTS.threshold_change);
            totalVolume=totalVolume+totalVolumeCurrent;
        end
    end
    
end


end