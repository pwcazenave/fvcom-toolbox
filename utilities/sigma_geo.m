function dist = sigma_geo(nlev, p_sigma)
% Generate a geometric sigma coordinate distribution.
%
% Mobj = sigma_gen(nlev, p_sigma)
%
% DESCRIPTION:
%   Generate a geometric vertical sigma coordinate distribution.
%
% INPUT:
%   nlev:       Number of sigma levels (layers + 1)
%   p_sigma:    1 for uniform sigma layers, 2 for parabolic function. See
%               page 308-309 in the FVCOM manual for examples.
%
% OUTPUT:
%   dist:       Geometric vertical sigma coordinate distribution.
%
% EXAMPLE USAGE:
%   Mobj = read_sigma(21, 2.0)
%
% Author(s):
%   Geoff Cowles (University of Massachusetts Dartmouth)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2013-04-23 Added help on the function and reformatted the code to
%   remove the FORTRAN in the else block.

dist = nan(1, nlev);

kb = nlev;

if p_sigma == 1
    for k = 1:nlev
        dist(k) = -((k-1)/(kb-1))^p_sigma;
    end
else
    for k = 1:(kb+1)/2
        dist(k) = -((k-1)/((kb+1)/2-1))^p_sigma/2;
    end
    for k = (kb+1)/2+1:kb
        dist(k) = ((kb-k)/((kb+1)/2-1))^p_sigma/2-1.0;
    end
end
