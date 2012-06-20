% wave friction using Soulsby page 80

function [tauw,fw] = wave_fric(u_orb,t_orb,z0,nu,rho)


ks = z0*30.; %Nikuradse roughness 

A = u_orb*t_orb/(2.*pi);   %semi-orbital excursion
r = A/ks; %relative roughness
Rw = u_orb*A/nu;  %wave Reynolds

fwr = 1.39*(A/z0)^(-.52); %rough turbulent flow friction fac
if(Rw > 5e5)
  B = .0521; N = .187; %smooth turbulent
else
  B = 2; N = 0.5;  % laminar
end;
fws = B*Rw^(-N);%smooth bed friction factor

fw = max(fwr,fws); %wave friction factor

tauw = 0.5*rho*fw*u_orb*u_orb; %wave shear stress
