function mu = SW_Viscosity(T,S)
% SW_Viscosity    Dynamic viscosity of seawater
%=========================================================================
% USAGE:  mu = SW_Viscosity(T,S)
%
% DESCRIPTION:
%   Dynamic viscosity of seawater at atmospheric pressure (0.1 MPa) using 
%   Eq. (22) given in [1] which best fit the data of [2], [3] and [4]. 
%   The pure water viscosity equation is a best fit to the data of [5]. 
%   Values at temperature higher than the normal boiling temperature 
%   are calculated at the saturation pressure.
%
% INPUT:  (all must have same dimensions)
%   T = temperature [degree C] (ITS-90)
%   S = salinity    [g/kg] (reference-composition salinity)
%
% OUTPUT:
%   mu = dynamic viscosity  [kg/m s]
%
% AUTHOR:  
%   Mostafa H. Sharqawy 12-18-2009, MIT (mhamed@mit.edu)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
% 
% VALIDITY: 0 < T < 180 C and 0 < S < 150 g/kg;
% 
% ACCURACY: 1.5%
% 
% REFERENCES:
%   [1] Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and Water Treatment, 2009
%   [2] B. M. Fabuss, A. Korosi, and D. F. Othmer, J., Chem. Eng. Data 14(2), 192, 1969.
%   [3] J. D. Isdale, C. M. Spence, and J. S. Tudhope, Desalination, 10(4), 319 - 328, 1972
%   [4] F. J. Millero, The Sea, Vol. 5, 3 � 80, John Wiley, New York, 1974
%   [5] IAPWS release on the viscosity of ordinary water substance 2008
%=========================================================================

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin ~=2
    error('SW_Viscosity.m: Must pass 2 parameters')
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
        if T(i,j)<0 | T(i,j)>180
            disp('Temperature is out of range for Viscosity function 10<T<180 C');
        end
        if S(i,j)<0 | S(i,j)>150
            disp('Salinity is out of range for Viscosity function 0<S<150 g/kg');
        end
        
        %------
        % BEGIN
        %------
        S(i,j)=S(i,j)/1000;
        a1 = 1.5700386464E-01;a2 = 6.4992620050E+01;a3 = -9.1296496657E+01;a4 = 4.2844324477E-05;
        mu_w(i,j) = a4 + 1/(a1*(T(i,j)+a2)^2+a3);
        a5 = 1.5409136040E+00;a6 = 1.9981117208E-02;a7 = -9.5203865864E-05;
        a8 = 7.9739318223E+00;a9 = -7.5614568881E-02;a10 = 4.7237011074E-04;
        A = a5 + a6 * T(i,j) + a7 * T(i,j)^2;
        B = a8 + a9 * T(i,j) + a10* T(i,j)^2;
        mu(i,j) = mu_w(i,j)*(1 + A*S(i,j) + B*S(i,j)^2);
    end
end
end