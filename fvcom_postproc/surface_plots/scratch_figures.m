close all

plot(allIdx,testChange)
hold on

plot(allIdx(negIdx),negIdx(negIdx==1)*0,'.')
plot(allIdx(negIdx),negIdx(negIdx==1)*0,'r.')

plot(allIdx(idxStart_dn),negIdx(idxStart_dn)*-0.001,'g.')
plot(allIdx(idxEnd_dn),negIdx(idxEnd_dn)*0.001,'k.')

% plot(allIdx(idxStart_tc),negIdx(idxStart_tc)*-0.002,'gx')
% plot(allIdx(idxEnd_tc),negIdx(idxEnd_tc)*0.002,'kx')

% ylim([-5e-3 5e-3])
axis([150 200 -5e-3 5e-3])

%%

for i=1:length(allInc)
    figure(1)
    clf
    plot(testTime,squeeze(phChange(allInc(i),1,1:end-1)))
    hold on
    for j=1:size(longOnesIdx{allInc(i)},1)
        plot(testTime(longOnesIdx{allInc(i)}(j,:)),zeros(2,1),'g-x')
    end
    text(min(testTime),0,num2str(allInc(i)))
    ylim([-0.0003 0.0003])
    pause(0.1)
end

%%

for i=1:length(allInc)
    figure(1)
    clf
    plot(testTime,squeeze(phChange(allInc(1),1,1:end-1)))
    hold on
    for j=1:size(longOnesIdx{1},1)
        plot(testTime(longOnesIdx{1}(j,:)),zeros(2,1),'g-x')
    end
    text(min(testTime),0,num2str(allInc(1)))
    ylim([-0.0003 0.0003])
    pause(0.1)
end

%%

figure(1)
clf
% plot(testTime,testData(testIdx,2:end))
plot(testTime,testData(2:end))
hold on
for j=1:size(longOnesIdx,1)
    plot(testTime(longOnesIdx(j,:)),zeros(2,1),'g-x')
end
