function totalVolume = get_runs(plotOPTS,data,cellVolume,lengthThreshold,changeThreshold)
% GET_RUNS Finds runs of continuous change in some value beyond some length
% and magnitude thresholds.
% 
% We're just doing lots of time series analyses here, so provide only a
% time series as DATA.

% Make sure we have a totalVolume value to output even if this time step
% doesn't match the threshold conditions specified.
totalVolume=0;

% Make an array of the time indices.
timeIdx=1:length(data);
% Calculate the differences ...
dataChange=diff(data);
% and find where they're negative.
negIdx=dataChange<0;

% Use the diff of the negative indices to find the start and end of each
% run.
dn=diff(negIdx);
dnIdxStart=timeIdx(dn>0)+1; % add one to move start along to correct index
dnIdxEnd=timeIdx(dn<0);

% Check if the start and end are identical (i.e. we have a spike).
if (numel(dnIdxStart)==1 || numel(dnIdxEnd)==1)
    if dnIdxStart==dnIdxEnd
%         warning('Single spike in time series, so carry on.')
        longOnesIdx(:,1:2)=nan(1,2);
        return
    end
end
% Check for no values
if sum(dnIdxStart)==0 || sum(dnIdxEnd)==0
%     warning('No appropriate values here, so carry on.')
    longOnesIdx(:,1:2)=nan(1,2);
    return
end
% Check the first index in dnIdxEnd is larger than the dnIdxStart. 
if dnIdxEnd(1)<dnIdxStart(1)
    % Strip it out and adjust the dnIdxStart index accordingly.
    dnIdxEnd=dnIdxEnd(2:end);
    dnIdxStart=dnIdxStart(1:end-1);
end

% If arrays are different lengths, lop off the last one from the longer
% array.
if length(dnIdxEnd)>length(dnIdxStart)
    timeIndices=(dnIdxEnd(1:end-1)-dnIdxStart)>=lengthThreshold;
elseif length(dnIdxEnd)<length(dnIdxStart)
    timeIndices=(dnIdxEnd-dnIdxStart(1:end-1))>=lengthThreshold;
elseif length(dnIdxEnd)==length(dnIdxStart)
    if dnIdxStart(1)==1
        timeIndices=(dnIdxEnd(1:end-1)-dnIdxStart(2:end))>=lengthThreshold;
    elseif dnIdxEnd(1)==1
        timeIndices=(dnIdxEnd(2:end)-dnIdxStart(1:end-1))>=lengthThreshold;
    else
        timeIndices=(dnIdxEnd-dnIdxStart)>=lengthThreshold;
    end
else
    return
end

if sum(timeIndices~=0)
    longOnesIdx(:,1)=dnIdxStart(timeIndices);
    longOnesIdx(:,2)=dnIdxEnd(timeIndices);
else
    return
end    

totalVolume=0;
for jj=1:size(longOnesIdx,1)
    % Average change
%     changeMetric=mean(data(longOnesIdx(jj,1):longOnesIdx(jj,2)));
    % Maximum change (i.e. closest to zero in our case)
    changeMetric=max(data(longOnesIdx(jj,1):longOnesIdx(jj,2)));
    if changeMetric<changeThreshold
        totalVolume=totalVolume+sum(cellVolume(longOnesIdx(jj,1):longOnesIdx(jj,2)));
    end
end

end

