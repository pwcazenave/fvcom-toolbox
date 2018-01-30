function dist = sigma_tanh(nlev,dl,du)
% Generate a tanh sigma coordinate distribution.
%
% Mobj = sigma_tanh(nlev, dl, du)
%
% DESCRIPTION:
%   Generate a tanh vertical sigma coordinate distribution.
%
% INPUT:
%   nlev:       Number of sigma levels (layers + 1)
%   dl:         The lower depth boundary from the bottom, down to which the
%               coordinates are parallel with uniform thickness.
%   du:         The upper depth boundary from the surface, up to which the
%               coordinates are parallel with uniform thickness.
%
% OUTPUT:
%   dist:       Tanh vertical sigma coordinate distribution.
%
% EXAMPLE USAGE:
%   Mobj = read_sigma(nlev, dl, du)
%
% Author(s):
%   Geoff Cowles (University of Massachusetts Dartmouth)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2013-04-23 Added help on the function and reformatted the code.

dist = zeros(1, nlev);

for k = 1:nlev-1
    x1 = dl+du;
    x1 = x1*(nlev-1-k)/(nlev-1);
    x1 = x1-dl;
    x1 = tanh(x1);
    x2 = tanh(dl);
    x3 = x2+tanh(du);
    dist(k+1) = (x1+x2)/x3-1.0;
end
