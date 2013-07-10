function [rDir,rMag,uRes,vRes]=do_residual(u,v,dt)
% DO_RESIDUAL Takes the u and v vectors of a model output and calculates
% the long-term direction and magnitude for that data.
% 
%   [RDIR,RMAG,URES,VRES]=DO_RESIDUAL(U,V,DT) takes the residual direction
%   (RDIR) and magnitude RMAG) of the data in U and V sampled at interval
%   DT. URES and UDIR are the summed U and V positions (the raw data for a
%   progresive vector diagram). Direction output is in degrees, vector
%   magnitude in units/s.
% 
% Pierre Cazenave PML 20/03/2012.
% 

% Loosely based on my original dfsuResidual.m and processResidual function
% for DHI's MIKE21 software, which in turn were based on Dave Lambkin's
% residual analysis scripts.
% 
% TODO: Make it possible to specify the average for all layers (i.e. NZ is
% all layers). 

% Let's do it...

toSecFactor=24*60*60;

nElements=size(u,1);
nLayers=size(u,2);
nTimeSteps=size(u,3);

% Some tidal assumptions. This will need to change in areas in which the
% diurnal tide dominates over the semidiurnal. 
tideCycle=(12+(25/60))/24;
tideWindow=ceil(tideCycle/dt);
tideDuration=(mean((dt*nTimeSteps)-tideCycle)-mean(tideCycle))*toSecFactor;

% Preallocate outputs.
uRes=zeros(nElements,nLayers,nTimeSteps);
vRes=zeros(nElements,nLayers,nTimeSteps);
uSum=nan(nElements,nTimeSteps,nLayers);
vSum=nan(nElements,nTimeSteps,nLayers);
uStart=nan(nElements,nLayers);
vStart=nan(nElements,nLayers);
uEnd=nan(nElements,nLayers);
vEnd=nan(nElements,nLayers);

for hh=1:nLayers
    uSum(:,:,hh)=cumsum(squeeze(u(:,hh,:)),2);
    vSum(:,:,hh)=cumsum(squeeze(v(:,hh,:)),2);
    for ii=1:nTimeSteps;
        uRes(:,hh,ii)=uRes(:,hh,ii)+(uSum(:,ii,hh).*(dt*toSecFactor));
        vRes(:,hh,ii)=vRes(:,hh,ii)+(vSum(:,ii,hh).*(dt*toSecFactor));
    end
    uStart(:,hh)=mean(squeeze(uRes(:,hh,1:tideWindow)),2);
    vStart(:,hh)=mean(squeeze(vRes(:,hh,1:tideWindow)),2);
    uEnd(:,hh)=mean(squeeze(uRes(:,hh,end-tideWindow:end)),2);
    vEnd(:,hh)=mean(squeeze(vRes(:,hh,end-tideWindow:end)),2);
end

uDiff=uEnd-uStart;
vDiff=vEnd-vStart;

% Calculate direction and magnitude.
rDir=atan2(uDiff,vDiff)*(180/pi); % in degrees.
rMag=sqrt(uDiff.^2+vDiff.^2)/tideDuration; % in units/s.