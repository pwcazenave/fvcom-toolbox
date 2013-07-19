function [fetch] = get_fetch(f,uwind,vwind,depth) 
%
% Determine fetch for given Cartesian wind speed or stress components 
%
% function get_fetch(uwind,vwind,f) 
%
% DESCRIPTION:
%   Display fetch relationship from fetch object 
%
% INPUT 
%   f     = fetch structure 
%   uwind = wind U10 or stress or other x-component
%   vwind = wind y-component
%   depth = [optional] depth at the station (default = uses bathymetry)
%   
%
% OUTPUT:
%   fetch in meters for that wind stress 
%
% EXAMPLE USAGE
%
%   fetch = get_fetch(myfetch,10.,0.,2.0)
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

[~,nZeta] = size(f.fetch);

%-------------------------------------------------
% find nearest points in theta/zeta space 
%-------------------------------------------------

% wind angle (-pi < wind angle < pi)
wangle = atan2(vwind,uwind);
[~,itheta] = min( abs(wangle-f.theta));

% zeta
if(exist('depth'));
  myzeta = f.zobs + depth;
else
  myzeta = 0.0;
end;
[~,izeta] = min( abs(myzeta-f.zeta));

% set fetch
fetch = f.fetch(itheta,izeta);

