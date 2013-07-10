function show_fetch(f) 
%
% Display fetch relationship from fetch object 
%
% function show_fetch(f) 
%
% DESCRIPTION:
%   Display fetch relationship from fetch object 
%
% INPUT 
%   f = fetch structure 
%
% OUTPUT:
%   plots 
%
% EXAMPLE USAGE
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

%-------------------------------------------------
% set dimensions 
%-------------------------------------------------

[nTheta,nZeta] = size(f.fetch);

%-------------------------------------------------
% option to plot 
%-------------------------------------------------
figure
scatter(f.xmesh,f.ymesh,5,f.hmesh); hold on;
plot(f.xobs,f.yobs,'b+'); hold on;
for i=1:nZeta
  plot(f.xobs+f.fetch(:,i)'.*cos(f.theta+pi),f.yobs+f.fetch(:,i)'.*sin(f.theta+pi),'k');
end;
figure
for i=1:nZeta
  plot((180/pi)*f.theta,f.fetch(:,i)/1000.,'k'); hold on;
end;
xlabel('Cartesian wind angle: e.g. pi/2 is south wind');
ylabel('fetch in km');
