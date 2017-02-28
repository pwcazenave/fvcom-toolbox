% Calculate the fraction of the rotor swept area is occupying each sigma layer
%
% Example Usage:
%
% sigma_frac = turbine_area_sigma(H, Ht, r, sigLay, plot_fig)
%
% Input Parameters:     H  - mean sea level (m)
%                       Ht - height of turbine hub above seabed (m)
%                       r  - turbine rotor radius (m)
%                       sigLay - number of sigma layers (not levels) in the model
%                       plot_fig - flag to plot a figure (optional)
%
% Rory O'Hara Murray, 19-Nov-2014
%
function sigma_frac = turbine_area_sigma(H, Ht, r, sigLay, plot_fig)

assert(nargin >= 4, 'Not enough arguments.');
assert(isnumeric(sigLay) && sigLay - fix(sigLay) < eps, 'sigLay (4th parameter) must be an integer number of sigma layers.');

if nargin<5
    plot_fig = false;
end

% make sure the water depth can accomodte the turbine at the desired hub
% height
if H < r + Ht
    warning(['Warning, a warter depth of ' num2str(H) ' m is insufficent, assuming water depth is ' num2str(r+Ht+0.01) ' m']);
    H = r + Ht + 0.01;
end


% Turbine and sigma layer parameters
elev  = 0; % water elevation above/below MSL - change this to see how the sigma layer occupation fraction changes with the tide
depth = H + elev;

dT = depth - Ht; % turbine hub depth

assert(dT>r, 'Turbine will stick out of water');

dLay = depth./sigLay;
zLev = [0:-dLay:-depth]';

% what sigma layer is the hub in?
drsl = zLev+dT; % depth of hub relative to each sigma level
hub_sigma = sum(drsl>=0);

%% draw sigma levels / layers
if plot_fig
    figure
    plot([-r r], zLev*[1 1])
    xlabel('Distance (m)')
    ylabel('Depth (m)')
    title([num2str(depth, '%2.0f') ' m water depth'])
    
    % draw rotor area
    a=0;
    b=-dT;
    ang = 0:pi/20:2*pi;
    x=r*cos(ang);
    y=r*sin(ang);
    hold on
    plot(a+x,b+y, a, b, 'o')
end
%% what fraction of the rotor area is in each sigma layer?

% loop trough all segments below the hub
dBot=-drsl(-drsl>=0);% the minimum of this array is hub height above a sigma level
numBot=sum(dBot<=r);
segmentsBotCum = [];
for ii=1:numBot
    phi = acos(dBot(ii)/r);
    sector = phi*r*r;
    triBot(ii) = r*sin(phi)*dBot(ii);
    segmentsBotCum(ii) = sector-triBot(ii);
end

% loop through all the segments above the hub
dTop=flip(drsl(drsl>=0));
numTop=sum(dTop<=r);
segmentsTopCum = [];
for ii=1:numTop
    phi = acos(dTop(ii)/r);
    sector = phi*r*r;
    triTop(ii) = r*sin(phi)*dTop(ii);
    segmentsTopCum(ii) = sector-triTop(ii);
end


segmentsTopCum2 = flip(segmentsTopCum);
segmentsTop = segmentsTopCum2-[0 segmentsTopCum2(1:end-1)];
segmentsBot = segmentsBotCum-[segmentsBotCum(2:end) 0];

% check that there are segments below/above hub, i.e. whether the rotors
% actually span multiple sigma layers
% sig_cent is area of the rotors in the layer the hub is in
if numTop>0 & numBot>0
    sigCent = pi*r*r-segmentsTopCum(1)-segmentsBotCum(1);   
elseif numTop==0 & numBot>0
    sigCent = pi*r*r-segmentsBotCum(1);
elseif numTop>0 & numBot==0
    sigCent = pi*r*r-segmentsTopCum(1);
elseif numTop==0 & numBot==0 % entire rotor is confined to one sigma layer
    sigCent = pi*r*r;   
end
    
% if sigCent is zero then the hub must exactely co-inside with a sigma
% level (unusual...) check it's larger than a very small area (or zero)
if sigCent>0
    % if hub co-insides exactely with a sigma level
    segments = [segmentsTop sigCent segmentsBot];
    segments_frac = segments./(pi*r*r);
    sig_span = hub_sigma + [-length(segmentsTop) length(segmentsBot)];
    numSigma = numTop+numBot+1; % total number of occupied sigma layers
else
    segments = [segmentsTop segmentsBot];
    segments_frac = segments./(pi*r*r);
    sig_span = hub_sigma + [-length(segmentsTop) length(segmentsBot)-1];
    numSigma = numTop+numBot; % total number of occupied sigma layers
end

% work out the fraction of the turbine area in each sigma layer
sigma_frac = zeros(1,sigLay);
sigma_frac(sig_span(1):sig_span(2)) = segments_frac;

total_frac = sum(sigma_frac);

