function ST_example
% Example usage of Sediment Toolbox
%
% function example
%
% DESCRIPTION:
%    Demonstrate ST Toolbox functionality 
%
% INPUT:
%
% OUTPUT:
%    
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
close all;
fprintf('       phi          class       d(mm)    Dstar    wset(mm/s)  taucr (Pa)  erate x1e-3(kg/(m^2-s))\n')
i = 0;
for phi=-8:11
	i = i + 1;
	phiclass    = ST_wentworth(phi);
	d(i)        = ST_phi2d(phi);
	Dstar(i)    = ST_Dstar(d(i));
	Wset(i)     = ST_wset(d(i));
	Taucr(i)    = ST_taucr(d(i));
	erate(i)    = ST_erate(d(i));
	fprintf('%10d %20s %8.4f %8.2f %9.4f %8.3f %8.3f\n',phi,phiclass,d(i)*1000,Dstar(i),Wset(i)*1000,Taucr(i),1000*erate(i))
end;

loglog(d*1000,Taucr)
title('critical shear stress')
xlabel('Grain diameter (mm)')
ylabel('critical shear Pa')
axis([.01,10,.01,10])