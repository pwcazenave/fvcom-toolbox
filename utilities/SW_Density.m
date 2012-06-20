function rho = SW_Density(T,S)
% SW_Density    Density of sea water
%=========================================================================
% USAGE:  rho = SW_Density(T,S)
%
% DESCRIPTION:
%   Density of seawater at atmospheric pressure (0.1 MPa) using Eq. (8)
%   given by [1] which best fit the data of [2] and [3]. The pure water
%   density equation is a best fit to the data of [4]. 
%   Values at temperature higher than the normal boiling temperature are
%   calculated at the saturation pressure.
%
% INPUT:  (all must have same dimensions)
%   T = temperature [degree C] (ITS-90)
%   S = salinity    [g/kg] (reference-composition salinity)
%
% OUTPUT:
%   rho = density   [kg/m^3]
%
% AUTHOR:  
%   Mostafa H. Sharqawy 12-18-2009, MIT (mhamed@mit.edu)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
% 
% VALIDITY: 0 < T < 180 C; 0 < S < 160 g/kg;
% 
% ACCURACY: 0.1%
% 
% REFERENCES:
%   [1] Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and Water Treatment, 2009.
%   [2] Isdale, and Morris, Desalination, 10(4), 329, 1972.
%   [3] Millero and Poisson, Deep-Sea Research, 28A (6), 625, 1981
%   [4]	IAPWS release on the Thermodynamic properties of ordinary water
%   substance, 1996.        
%   UPDATED 09-23-2010 modified to now handle matrices and commented out
%   range checking.
%=========================================================================

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin ~=2
    error('SW_Density.m: Must pass 2 parameters')
end

% CHECK S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);

% CHECK THAT S & T HAVE SAME SHAPE
if (ms~=mt) | (ns~=nt)
    error('check_stp: S & T must have same dimensions')
end

% CHECK THAT S & T ARE WITHIN THE FUNCTION RANGE
vectorsize=size(S);
for i = 1:vectorsize(1,1)
    for j = 1:vectorsize(1,2)
%         if T(i,j)>180 | T(i,j)<0
%             disp('Temperature is out of range for density function 0 < T < 180 C');
%         end
%         if S(i,j)<0 | S(i,j)>160
%             disp('Salinity is out of range for density function 0 < S < 160 g/kg');
%         end
        
        %------
        % BEGIN
        %------
        s(i,j)=S(i,j)/1000;
        a1=9.9992293295E+02;a2=2.0341179217E-02;a3=-6.1624591598E-03;a4=2.2614664708E-05;a5=-4.6570659168E-08;
        b1=8.0200240891E+02;b2=-2.0005183488E+00;b3=1.6771024982E-02;b4=-3.0600536746E-05;b5=-1.6132224742E-05;
        rho_w(i,j) = a1 + a2*T(i,j) + a3*T(i,j)^2 + a4*T(i,j)^3 + a5*T(i,j)^4;
        D_rho(i,j) = b1*s(i,j) + b2*s(i,j)*T(i,j) + b3*s(i,j)*T(i,j)^2 + b4*s(i,j)*T(i,j)^3 + b5*s(i,j)^2*T(i,j)^2;
        rho(i,j) = rho_w(i,j) + D_rho(i,j);
    end
end;
end