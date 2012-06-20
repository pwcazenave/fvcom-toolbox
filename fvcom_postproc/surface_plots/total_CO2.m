% Get the total CO2 for a day of model.
close all

leakidx=1316;
dt=3600;

[xt,zt,ttot]=size(FVCOM.DYE);

nt=24; % Let's do a day.
nBuff=12; % half day buffer
dtStartIdx=(3600*nBuff)/dt; % in case we don't have hourly sampled output
dtEndIdx=((3600*nt)/dt)+dtStartIdx;

% 3D array of the layer thicknesses
thickFactor=repmat(abs(diff(FVCOM.siglev,1,2)),[1,1,ttot]);
% 3D array of the total depths at each element (tide + still water depth)
totalDepth=permute(repmat(FVCOM.zeta+repmat(FVCOM.h,1,ttot),[1,1,zt]),[1,3,2]);
% Now we have two 3D arrays, we just need to multiply totalDepth by
% thickFactor to get layer thickness in metres
layerThickness=thickFactor.*totalDepth;

% Now, integrate the CO2 by depth across a day's worth of time steps.
