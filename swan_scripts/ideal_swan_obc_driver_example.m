function ideal_swan_obc_driver();
	
% Example Driver to Generate idealized SWAN boundary forcing using ideal_swan_obc
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

t1 = greg2mjulian(2007,1,1,0,0,0);
t2 = greg2mjulian(2007,2,1,0,0,0);
time = t1:(6./24):t2;
hs = 2.*sin((time-time(1))*2*pi/15);
tp = 15. + 5.*sin((time-time(1))*2*pi/15);
dir = 180*ones(prod(size(time)),1);
spread = 10.;
nodefile = 'tst.node';
inc = 5;
ideal_swan_obc(nodefile,time,hs,tp,dir,inc,spread);  

