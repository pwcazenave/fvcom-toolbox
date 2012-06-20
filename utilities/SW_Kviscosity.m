function new = SW_Kviscosity(T,S)
% Calculates the kinemtic viscosity [m^2/s] from the dynamic viscosity and
% density functions
mu = SW_Viscosity(T,S);
rho = SW_Density(T,S);
new = mu./rho;
end