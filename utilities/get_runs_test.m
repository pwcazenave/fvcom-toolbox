% Figure out a way to identify runs longer than some specified length for
% the purposes of calculating the volume of water which is subjected to a
% change in a parameter of some value. 

%% 1D

% close all

clear testChange negIdx dn startEndIdx gotTo inc longOnesIdx

testIdx=88;

dtJump=plotOPTS.change_type;
lengthThreshold=12;
changeThreshold=-0.5;

% Model results
testData=(squeeze(phChange(testIdx,1,:)));
testTime=plotOPTS.Time_record(1:end-dtJump-1);
% Synthetic data
% testData=[10,10,9,8,7,6,4,2,0,-2,-4,-6,-8,-10,-20,-30,-40,-50,-60,-70,-80,-90,-100,-110,-109,-108,-107,-106,-105,-104,-103,-102,-101,-100,-99,-97,-95,-96,-97,-98,-99,-99,-99,-99,-99,-100,-101,-102,-103,-105,-107,-109,-111,-113,-115,-117,-119,-129,-130,-130,-129,-128,-127];
% testTime=1:length(testData)-1;

testChange=diff(testData);

allIdx=1:length(testData); % array of index positions
negIdx=testChange<0; % where the negative values are. 
dn=diff(negIdx); % the boundaries of those negative data

% Find locations at which dn changes from zero to non-zero
dnIdxStart=allIdx(dn==1)+1; % add one to move start along to correct index
dnIdxEnd=allIdx(dn==-1);

% Check if the start and end are identical (i.e. we have a spike).
if (numel(dnIdxStart)==1 || numel(dnIdxEnd)==1)
    if dnIdxStart==dnIdxEnd
        warning('Single spike in time series, so carry on.')
        continue
    end
end

if sum(dnIdxStart)==0 || sum(dnIdxEnd)==0
    % No appropriate values here, so carry on.
    continue
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

% Check the change magnitude for those time steps. Here we have to decide
% if we're doing average change over the n-hour time period, or the maximum
% instantaneous change, or ... For a first go, we'll do the average change
% and use that for the threshold comparison.
checkInc=1;
totalVolume=0;
for jj=1:size(longOnesIdx,2)
    averageChange=mean(testData(longOnesIdx(jj,1):longOnesIdx(jj,2)));
    if averageChange<changeThreshold
        totalVolume=totalVolume+cellVolume(longOnesIdx(jj,1):longOnesIdx(jj,2));
        checkInc=checkInc+1;
    end
end

close all

figure(2)
plot(testTime,testChange)
hold on
plot(testTime(testChange<0),testChange(testChange<0),'r.')
% plot(testTime(4:end)-((testTime(2)-testTime(1))/2),dn*0.01,'gp')
for i=1:size(longOnesIdx,1)
    plot(testTime(longOnesIdx(i,:)),zeros(1,length(longOnesIdx(i,:))),'g-x')
    text(mean(testTime(longOnesIdx(i,:))),0.01,num2str((testTime(longOnesIdx(i,2))-testTime(longOnesIdx(i,1)))*24))
end

%% 1Dv2

close all

clear testChange negIdx dn startEndIdx gotTo inc longOnesIdx

testIdx=124;

dtJump=plotOPTS.change_type;
lengthThreshold=4;
changeThreshold=-1.5;

% Model results
testData=(squeeze(phChange(testIdx,1,:)));
testTime=plotOPTS.Time_record(1:end-dtJump-1);
% Synthetic data
% testData=[10,10,9,8,7,6,4,2,0,-2,-4,-6,-8,-10,-20,-30,-40,-50,-60,-70,-80,-90,-100,-110,-109,-108,-107,-106,-105,-104,-103,-102,-101,-100,-99,-97,-95,-96,-97,-98,-99,-99,-99,-99,-99,-100,-101,-102,-103,-105,-107,-109,-111,-113,-115,-117,-119,-129,-130,-130,-129,-128,-127];
% testTime=1:length(testData)-1;

testChange=diff(testData);

allIdx=1:length(testData); % array of index positions
negIdx=testChange<0; % where the negative values are. 
dn=diff(negIdx); % the boundaries of those negative data

% Find locations at which dn changes from zero to non-zero
dnIdxStart=allIdx(dn==1)+1; % add one to move start along to correct index
dnIdxEnd=allIdx(dn==-1);

% Check if the start and end are identical (i.e. we have a spike).
if (numel(dnIdxStart)==1 || numel(dnIdxEnd)==1)
    if dnIdxStart==dnIdxEnd
        warning('Single spike in time series, so carry on.')
    end
end
% Check for no values
if sum(dnIdxStart)==0 || sum(dnIdxEnd)==0
    % No appropriate values here, so carry on.
    warning('No appropriate values here, so carry on.')
end
% Check the first index in dnIdxEnd is larger than the dnIdxStart. 
if dnIdxEnd(1)<dnIdxStart(1)
    % Strip it out and adjust the dnIdxStart index accordingly.
    dnIdxEnd=dnIdxEnd(2:end);
    dnIdxStart=dnIdxStart(1:end-1);
end


% If arrays are different lengths, lop off the last one from the longer
% array.
if size(dnIdxEnd,2)>size(dnIdxStart,2)
    longOnesIdx(:,1)=dnIdxStart((dnIdxEnd(1:end-1)-dnIdxStart)>=lengthThreshold);
    longOnesIdx(:,2)=dnIdxEnd((dnIdxEnd(1:end-1)-dnIdxStart)>=lengthThreshold);
elseif size(dnIdxEnd,2)<size(dnIdxStart,2)
    longOnesIdx(:,1)=dnIdxStart((dnIdxEnd-dnIdxStart(1:end-1))>=lengthThreshold);
    longOnesIdx(:,2)=dnIdxEnd((dnIdxEnd-dnIdxStart(1:end-1))>=lengthThreshold);
else
    if dnIdxStart(1)==1
        longOnesIdx(:,1)=dnIdxStart((dnIdxEnd(1:end-1)-dnIdxStart(2:end))>=lengthThreshold);
        longOnesIdx(:,2)=dnIdxEnd((dnIdxEnd(1:end-1)-dnIdxStart(2:end))>=lengthThreshold);
    elseif dnIdxEnd(1)==1
        longOnesIdx(:,1)=dnIdxStart((dnIdxEnd(2:end)-dnIdxStart(1:end-1))>=lengthThreshold);
        longOnesIdx(:,2)=dnIdxEnd((dnIdxEnd(2:end)-dnIdxStart(1:end-1))>=lengthThreshold);
    else
        longOnesIdx(:,1)=dnIdxStart((dnIdxEnd-dnIdxStart)>=lengthThreshold);
        longOnesIdx(:,2)=dnIdxEnd((dnIdxEnd-dnIdxStart)>=lengthThreshold);
    end        
end

close all

plot(testTime,testChange)
hold on
plot(testTime(testChange<0),testChange(testChange<0),'r.')
% plot(testTime(4:end)-((testTime(2)-testTime(1))/2),dn*0.01,'gp')
for i=1:size(longOnesIdx,1)
    plot(testTime(longOnesIdx(i,:)),zeros(1,length(longOnesIdx(i,:))),'g-x')
    text(mean(testTime(longOnesIdx(i,:))),0.01,num2str((testTime(longOnesIdx(i,2))-testTime(longOnesIdx(i,1)))*24))
end


%% 2D

clear testChange negIdx dn startEndIdx gotTo inc allInc longOnesIdx

dtJump=plotOPTS.change_type;
lengthThreshold=7;
changeThreshold=-0.1;

testData=squeeze(phChange(:,1,:));
testTime=plotOPTS.Time_record(1:end-dtJump-1);

inc=1;

% Easiest way I can think of doing this is to just loop through each
% element, doing each timestep as a 1D analysis.
for i=1:size(testData,1);
    disp(i)
    testChange=diff(testData(i,:),1,2);
    allIdx=1:length(testChange);
    negIdx=testChange<0; % where the negative values are. 
    dn=diff(negIdx,1,2);

    % Add one to move start along to correct index. Since the diff
    % describes the relationship between the first and second value in
    % testData at the first value of dn. So, if you have a change between
    % the first and second testData values, then the start of the run is
    % necessarily at the second index in testData. For the end though, the
    % index in dn will describe the change from the run to a non-run part,
    % so we don't need to adjust that in the same way. We need also to
    % adjust the values of testChange when comparine them with 
    %
    % allIdx:   |  1  |  2  |  3  |  4  |  5  |  6  |  7  |
    % testData: |  P  |     |     |     |     |  P  |     |
    %           |     |  N  |  N  |  N  |  N  |     |  N  | 
    % negIdx:   |  0  |  1  |  1  |  1  |  1  |  0  |  1  | 
    % dn:       |  1  |  0  |  0  |  0  | -1  |  1  |  ~  |
    % 
    % So you can see that dn is positive before the start of the N run
    % whilst dn is negative at the last element of the run. So we need to
    % increase the idxStart value by one. The idxEnd doesn't need to be
    % modified this way.
    idxStart_dn=allIdx(dn==1)+1;
%     idxStart_tc=allIdx(testChange(1:end-1)<changeThreshold);
%     dnIdxStart=intersect(idxStart_dn,idxStart_tc);
    dnIdxStart=idxStart_dn;
    idxEnd_dn=allIdx(dn==-1);
%     idxEnd_tc=allIdx(testChange(1:end-1)<changeThreshold);
%     dnIdxEnd=intersect(idxEnd_dn,idxEnd_tc);
    dnIdxEnd=idxEnd_dn;

    % QC the picked results.

    % Check if the start and end are identical (i.e. we have a spike).
    if (numel(dnIdxStart)==1 || numel(dnIdxEnd)==1)
        if dnIdxStart==dnIdxEnd
            warning('Single spike in time series, so carry on.')
            continue
        end
    end
    % Check for no values
    if sum(dnIdxStart)==0 || sum(dnIdxEnd)==0
        % No appropriate values here, so carry on.
        warning('No appropriate values here, so carry on.')
        continue
    end
    % Check the first index in dnIdxEnd is larger than the dnIdxStart. 
    if dnIdxEnd(1)<dnIdxStart(1)
        % Strip it out and adjust the dnIdxStart index accordingly.
        dnIdxEnd=dnIdxEnd(2:end);
        dnIdxStart=dnIdxStart(1:end-1);
    end

    % If arrays are different lengths, lop off the last one from the longer
    % array.
    if size(dnIdxEnd,2)>size(dnIdxStart,2)
        % Not ideal with the try's, but it seems to work...
        try
            timeRange(:,1)=dnIdxStart((dnIdxEnd(1:end-1)-dnIdxStart)>=lengthThreshold);
            timeRange(:,2)=dnIdxEnd((dnIdxEnd(1:end-1)-dnIdxStart)>=lengthThreshold);
        catch %#ok<CTCH>
            continue
        end
    elseif size(dnIdxEnd,2)<size(dnIdxStart,2)
        try
            timeRange(:,1)=dnIdxStart((dnIdxEnd-dnIdxStart(1:end-1))>=lengthThreshold);
            timeRange(:,2)=dnIdxEnd((dnIdxEnd-dnIdxStart(1:end-1))>=lengthThreshold);
        catch %#ok<CTCH>
            continue
        end
    else
        try
            timeRange(:,1)=dnIdxStart((dnIdxEnd-dnIdxStart)>=lengthThreshold);
            timeRange(:,2)=dnIdxEnd((dnIdxEnd-dnIdxStart)>=lengthThreshold);
        catch %#ok<CTCH>
            continue
        end
    end
    
    % Now do the magnitude check now.
    timeRange=timeRange(testChange(timeRange(:,1))<changeThreshold & testChange(timeRange(:,2))<changeThreshold,:);
    % If the start and end ranges are non-zero in length, add them to the
    % longOnesIdx cell array.
    if numel(timeRange)>=2
        longOnesIdx{inc}=timeRange;
        allInc(inc)=i;
        inc=inc+1;
    end
end

%% 2Dv2

clear testChange negIdx dn startEndIdx gotTo inc allInc longOnesIdx timeRange

dtJump=plotOPTS.change_type;
lengthThreshold=2;
changeThreshold=-1;

testData=squeeze(phChange(:,1,:));
testTime=plotOPTS.Time_record(1:end-dtJump-1);

inc=1;

% Easiest way I can think of doing this is to just loop through each
% element, doing each timestep as a 1D analysis.
for i=222%1:size(testData,1);
    disp(i)
    % Model results
    testData=squeeze(phChange(i,1,:))*10000000000;
    testTime=plotOPTS.Time_record(1:end-dtJump-1);
    % Synthetic data
    % testData=[10,10,9,8,7,6,4,2,0,-2,-4,-6,-8,-10,-20,-30,-40,-50,-60,-70,-80,-90,-100,-110,-109,-108,-107,-106,-105,-104,-103,-102,-101,-100,-99,-97,-95,-96,-97,-98,-99,-99,-99,-99,-99,-100,-101,-102,-103,-105,-107,-109,-111,-113,-115,-117,-119,-129,-130,-130,-129,-128,-127];
    % testTime=1:length(testData)-1;

    testChange=diff(testData);

    allIdx=1:length(testData); % array of index positions
    negIdx=testChange<0; % where the negative values are. 
    dn=diff(negIdx); % the boundaries of those negative data

    % Find indices at which dn changes from zero to non-zero
    dnIdxStart=allIdx(dn==1)+1; % add one to move start along to correct index
    dnIdxEnd=allIdx(dn==-1);

    % Check if the start and end are identical (i.e. we have a spike).
    if (numel(dnIdxStart)==1 || numel(dnIdxEnd)==1)
        if dnIdxStart==dnIdxEnd
            warning('Single spike in time series, so carry on.')
            continue
        end
    end
    % Check for no values
    if sum(dnIdxStart)==0 || sum(dnIdxEnd)==0
        % No appropriate values here, so carry on.
        warning('No appropriate values here, so carry on.')
        continue
    end
    % Check the first index in dnIdxEnd is larger than the dnIdxStart. 
    if dnIdxEnd(1)<dnIdxStart(1)
        % Strip it out and adjust the dnIdxStart index accordingly.
        dnIdxEnd=dnIdxEnd(2:end);
        dnIdxStart=dnIdxStart(1:end-1);
    end

    % If arrays are different lengths, lop off the last one from the longer
    % array.
    clear timeRange
    if size(dnIdxEnd,2)>size(dnIdxStart,2)
        % Not ideal with the try's, but it seems to work...
        try
            timeRange(:,1)=dnIdxStart((dnIdxEnd(1:end-1)-dnIdxStart)>=lengthThreshold);
            timeRange(:,2)=dnIdxEnd((dnIdxEnd(1:end-1)-dnIdxStart)>=lengthThreshold);
        catch %#ok<CTCH>
            continue
        end
    elseif size(dnIdxEnd,2)<size(dnIdxStart,2)
        try
            timeRange(:,1)=dnIdxStart((dnIdxEnd-dnIdxStart(1:end-1))>=lengthThreshold);
            timeRange(:,2)=dnIdxEnd((dnIdxEnd-dnIdxStart(1:end-1))>=lengthThreshold);
        catch %#ok<CTCH>
            continue
        end
    else
        try
            timeRange(:,1)=dnIdxStart((dnIdxEnd-dnIdxStart)>=lengthThreshold);
            timeRange(:,2)=dnIdxEnd((dnIdxEnd-dnIdxStart)>=lengthThreshold);
        catch %#ok<CTCH>
            continue
        end
    end

    % Now do the magnitude check now.
    timeRange=timeRange(testChange(timeRange(:,1))<changeThreshold & testChange(timeRange(:,2))<changeThreshold,:);
    
    % If the start and end ranges are non-zero in length, add them to the
    % longOnesIdx cell array.
    if numel(timeRange)>=2
        longOnesIdx{inc}=timeRange;
        allInc(inc)=i;
        inc=inc+1;
    end

%     if size(dnIdxEnd,2)>size(dnIdxStart,2)
%         try
%             timeRange
%             longOnesIdx{inc}(:,1)=dnIdxStart((dnIdxEnd(1:end-1)-dnIdxStart)>=lengthThreshold);
%             longOnesIdx{inc}(:,2)=dnIdxEnd((dnIdxEnd(1:end-1)-dnIdxStart)>=lengthThreshold);
%             allInc(inc)=i;
%             inc=inc+1;
%         end
%     elseif size(dnIdxEnd,2)<size(dnIdxStart,2)
%         try
%             longOnesIdx{inc}(:,1)=dnIdxStart((dnIdxEnd-dnIdxStart(1:end-1))>=lengthThreshold);
%             longOnesIdx{inc}(:,2)=dnIdxEnd((dnIdxEnd-dnIdxStart(1:end-1))>=lengthThreshold);
%             allInc(inc)=i;
%             inc=inc+1;
%         end
%     else
%         try
%             longOnesIdx{inc}(:,1)=dnIdxStart((dnIdxEnd-dnIdxStart)>=lengthThreshold);
%             longOnesIdx{inc}(:,2)=dnIdxEnd((dnIdxEnd-dnIdxStart)>=lengthThreshold);
%             allInc(inc)=i;
%             inc=inc+1;
%         end
%     end

   
%     close all
% 
%     plot(testTime,testChange)
%     hold on
%     plot(testTime(testChange<0),testChange(testChange<0),'r.')
%     % plot(testTime(4:end)-((testTime(2)-testTime(1))/2),dn*0.01,'gp')
%     for j=1:size(longOnesIdx{inc-1},1)
%         plot(testTime(longOnesIdx{inc}(j,:)),zeros(1,length(longOnesIdx{inc}(j,:))),'g-x')
%         text(mean(testTime(longOnesIdx{inc-1}(j,:))),0.01,num2str((testTime(longOnesIdx{inc-1}(j,2))-testTime(longOnesIdx{inc-1}(j,1)))*24))
%     end

end
%%
close all

[tx,ty]=meshgrid(1:size(testChange,2),1:size(testChange,1));

pcolor(tx,ty,testChange); shading flat; axis tight
colorbar
caxis([-1 1])
hold on
plot(tx(allIdx(dn<0)),ty(allIdx(dn<0)),'r.')
plot(tx(allIdx(dn>0)),ty(allIdx(dn>0)),'g.')
for i=1:size(longOnesIdx,1)
    plot(tx(longOnesIdx(i,1)),ty(longOnesIdx(i,2)),'g-x')
end

%% 1D (old way)

% % Find continuous negatives for more than changeThreshold "hours".
% inc=0;
% gotTo=1;
% for i=1:length(dn)
%     % Check whether we're in a flat big (i.e. the change is constant) and
%     % whether the value of the change exceeds the threshold change
%     % specified.
%     if dn(i)==0 && testChange(i)<0
%         fprintf('dn: %i\ttestChange: %i\tinc: %i\n',dn(i),testChange(i),inc)
%         % Check what the value of inc is:
%         % - If it's zero, we're at the start of a run.
%         % - If it's non-zero, we're part way through a run.
%         % - If it's negative, the last value was the end of a run.
%         if inc==0
%             % Start of a run
%             startEndIdx(gotTo,1)=i;
%             gotTo=gotTo+1;
%             inc=1;
%         elseif inc>0
%             % Midway through a run
%             inc=inc+1;
%         else
%             warning('Uhhh... this shouldn''t ever appear...')
%         end
%     elseif dn(i)~=0 && testChange(i)<0
%         % At the end of a run when diff is non-zero and increment is
%         % greater than zero.
%         if inc>0
%             startEndIdx(gotTo-1,2)=i;
%             inc=0;
%         end
%     end
%     
% end


