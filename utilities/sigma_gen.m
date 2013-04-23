function dist = sigma_gen(nlev,dl,du,kl,ku,zkl,zku,h,hmin)
% Generate a generalised sigma coordinate distribution.
%
% Mobj = sigma_gen(nlev, dl, du, kl, ku, zkl, zku, h, hmin)
%
% DESCRIPTION:
%   Generate a uniform or hybrid vertical sigma coordinate system.
%
% INPUT:
%   nlev:       Number of sigma levels (layers + 1)
%   dl:         The lower depth boundary from the bottom, down to which the
%               coordinates are parallel with uniform thickness.
%   du:         The upper depth boundary from the surface, up to which the
%               coordinates are parallel with uniform thickness.
%   kl:         ?
%   ku:         ?
%   zkl:        ?
%   zku:        ?
%   h:          Water depth.
%   hmin:       Minimum water depth.
%
% OUTPUT:
%   dist:       Generalised vertical sigma coordinate distribution.
%
% EXAMPLE USAGE:
%   Mobj = sigma_gen(nlev, dl, du, kl, ku, zkl, zku, h, hmin)
%
% Author(s):
%   Geoff Cowles (University of Massachusetts Dartmouth)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2013-04-23 Added help on the function and reformatted the code.

dist = nan(1, nlev);

if h < hmin
    dist(1) = 0.0;
    dl2 = 0.001;
    du2 = 0.001;
    for k = 1:nlev-1
        x1 = dl2+du2;
        x1 = x1*double(nlev-1-k)/double(nlev-1);
        x1 = x1 - dl2;
        x1 = tanh(x1);
        x2 = tanh(dl2);
        x3 = x2+tanh(du2);
        dist(k+1) = (x1+x2)/x3-1.0;
    end
else
    %dr=(h-sum(zku)-sum(zkl))/h/double(nlev-ku-kl-1);
    dr = (h-du-dl)/h/double(nlev-ku-kl-1);
    dist(1) = 0.0;

    for k = 2:ku+1
        dist(k) = dist(k-1)-zku(k-1)/h;
        %     fprintf('building z %f %f %f %f \n',z(k),zku(k-1),h,zku(k-1)/h)
    end

    for k = ku+2:nlev-kl
        dist(k) = dist(k-1)-dr;
        %     fprintf('building z %f %f \n',z(k),dr)
    end

    kk = 0;
    for k = nlev-kl+1:nlev
        kk = kk+1;
        dist(k) = dist(k-1)-zkl(kk)/h;
        %     fprintf('building z %f %f \n',z(k),zkl(kk))
    end
end
