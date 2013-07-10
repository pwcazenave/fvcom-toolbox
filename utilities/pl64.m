function xf=pl66tn(x,dt,T);
% low pass filter (33 hr) 
% PL66TN: pl66t for variable dt and T
% xf=PL66TN(x,dt,T) computes low-passed series xf from x
% using pl66 filter, with optional sample interval dt (hrs)
% and filter half amplitude period T (hrs) as input for
% non-hourly series.
%
% INPUT:  x=time series (must be column array)
%         dt=sample interval time [hrs] (Default dt=1)
%         T=filter half-amp period [hrs] (Default T=33)
%
% OUTPUT: xf=filtered series

% NOTE: both pl64 and pl66 have the same 33 hr filter
% half-amplitude period. pl66 includes additional filter weights
% upto and including the fourth zero crossing at 2*T hrs.

% The PL64 filter is described on p. 21, Rosenfeld (1983), WHOI
% Technical Report 85-35. Filter half amplitude period = 33 hrs.,
% half power period = 38 hrs. The time series x is folded over
% and cosine tapered at each end to return a filtered time series
% xf of the same length. Assumes length of x greater than 132.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10/30/00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default to pl64
if (nargin==1); dt=1; T=33; end

cutoff=T/dt;
fq=1./cutoff;
nw=2*T./dt;
nw=round(nw);
%disp(['number of weights = ',int2str(nw)])
nw2=2.*nw;

[npts,ncol]=size(x);
if (npts<ncol);x=x';[npts,ncol]=size(x);end
xf=x;

% generate filter weights
j=1:nw;
t=pi.*j;
den=fq.*fq.*t.^3;
wts=(2.*sin(2.*fq.*t)-sin(fq.*t)-sin(3.*fq.*t))./den;
% make symmetric filter weights
wts=[wts(nw:-1:1),2.*fq,wts];
wts=wts./sum(wts);% normalize to exactly one
% plot(wts);grid;
% title(['pl64t filter weights for dt = ',num2str(dt),' and T = ',num2str(T)])
% xlabel(['number of weights = ',int2str(nw)]);pause;

% fold tapered time series on each end
cs=cos(t'./nw2);
jm=[nw:-1:1];

for ic=1:ncol
% ['column #',num2str(ic)]
 jgd=find(~isnan(x(:,ic)));
 npts=length(jgd);
 if (npts>nw2)
%detrend time series, then add trend back after filtering
  xdt=detrend(x(jgd,ic));
  trnd=x(jgd,ic)-xdt;
  y=[cs(jm).*xdt(jm);xdt;cs(j).*xdt(npts+1-j)];
% filter
  yf=filter(wts,1.0,y);
% strip off extra points
  xf(jgd,ic)=yf(nw2+1:npts+nw2);
% add back trend
  xf(jgd,ic)=xf(jgd,ic)+trnd;
 else
 'warning time series is too short'
 end
end
