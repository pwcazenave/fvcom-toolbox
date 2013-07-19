% Daily CO2 input for bottom cell

close all
clc

leakidx=1316;
baseidx=1;
dt=3600;

% Get bottom element input for 24 hours (skipping the buffer at the
% beginning of half a day or so)
nt=24; % Let's do a day.
nBuff=12; % half day buffer
dtStartIdx=(3600*nBuff)/dt; % in case we don't have hourly sampled output
dtEndIdx=((3600*nt)/dt)+dtStartIdx;

% Should just be the input amount at the bottom cell times the time step
% times the number of time steps in a day.
fprintf('Daily input (pica): %g\n',sum(squeeze(FVCOM.DYE(leakidx,baseidx,dtStartIdx:dtEndIdx))*dt))


% My version: ((c*dt*n)/v)*86400
fprintf('Daily input (pica): %g\n',(sum(squeeze(FVCOM.DYE(leakidx,baseidx,dtStartIdx:dtEndIdx)))*dt/(FVCOM.art1(leakidx)*abs(diff(FVCOM.siglev(leakidx,baseidx:baseidx+1)))))*86400)
% Riqui's version: c*dt*n*v
fprintf('Daily input (rito): %g\n',sum(squeeze(FVCOM.DYE(leakidx,baseidx,dtStartIdx:dtEndIdx)))*dt*FVCOM.art1(leakidx)*abs(diff(FVCOM.siglev(leakidx,baseidx:baseidx+1))))

%% A new dawn...

home

% Modelled
modelledInput=zeros(dtEndIdx-dtStartIdx+1,1);
dtIter=dtStartIdx:dtEndIdx;
for i=1:length(dtIter)
    cellVolume=FVCOM.art1(leakidx)*((FVCOM.h(leakidx)+FVCOM.zeta(leakidx,dtIter(i)))*abs(diff(FVCOM.siglev(leakidx,baseidx:baseidx+1))));
    if i>1
        modelledInput(i)=modelledInput(i-1)+(((FVCOM.DYE(leakidx,baseidx,dtIter(i))*dt)/cellVolume));
    end
end
TCO2=cumsum((squeeze(FVCOM.DYE(leakidx,baseidx,dtStartIdx:dtEndIdx))*dt)/cellVolume);

% Theoretical
theoreticalInput=zeros(dtEndIdx-dtStartIdx+1,1);
for i=1:length(dtIter)
    % Use the real cell volumes
    cellVolume=FVCOM.art1(leakidx)*((FVCOM.h(leakidx)+FVCOM.zeta(leakidx,dtIter(i)))*abs(diff(FVCOM.siglev(leakidx,baseidx:baseidx+1))));
    if i>1
        theoreticalInput(i)=theoreticalInput(i-1)+(((500*dt)/cellVolume));
    end
end

fprintf('Daily input (pica2): %g\n',modelledInput(end))

fprintf('Daily input (theory): %g\n',theoreticalInput(end))
fprintf('Daily input (theory2): %g\n',TCO2(end))


%%

% Compare the theoretical input (i.e. the input rate by the number of time
% steps) against the modelled input at the leak point.
close all

clear modelledTCO2 theoreticalTCO2

leakidx=1316;
baseidx=1;
dt=3600;
inputCO2=500;

colours=hsv(length(baseidx));

nt=24; % Let's do a day.
nBuff=12; % half day buffer
dtStartIdx=(3600*nBuff)/dt; % in case we don't have hourly sampled output
dtEndIdx=((3600*nt)/dt)+dtStartIdx;

% Get the element area
elemArea=FVCOM.art1(leakidx);

% Get the modelled values
for i=1:length(baseidx);
    % Get the height of the base element
    elemHeight=(FVCOM.zeta(leakidx,dtStartIdx:dtEndIdx)+FVCOM.h(leakidx)).*abs(diff(FVCOM.siglev(leakidx,baseidx(i):baseidx(i)+1)));
    % Get the modelled CO2, and adjust for volume of base element
    modelledTCO2(:,i)=(cumsum(squeeze(FVCOM.DYE(leakidx,baseidx(i),dtStartIdx:dtEndIdx)))./(elemArea.*elemHeight'))*dt;
end

% And the theoretical values (with the modelled tide effect i.e. change in
% volume).
theoreticalTCO2=((1:dtEndIdx-dtStartIdx+1)*500./(elemArea.*elemHeight))*dt;

figure(1)
hold on
for i=1:length(baseidx);
    plot(Time_record(dtStartIdx:dtEndIdx),modelledTCO2,'-x','LineWidth',3)
end
plot(Time_record(dtStartIdx:dtEndIdx),theoreticalTCO2,'r--x','LineWidth',2)
axis('tight')
ylabel('DYE')
xlabel('Time (days)')
legend('Modelled','Theoretical','Location','NorthWest')
legend('BoxOff')