function [Mobj] = sigma_tanh(nlev,dl,du,Mobj,sigma_file)
% Generate a tanh sigma coordinate distribution.
%
% Mobj = sigma_tanh(nlev, dl, du,Mobj)
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
%   Mobj:       [optional] mesh object file
%
% OUTPUT:
%   Mobj.sigma:       Tanh vertical sigma coordinate distribution.
%   Mobj.sig...     All the sigma layers variables such as siglev, siglevz,
%                   siglayz, etc...
%
% EXAMPLE USAGE:
%   Mobj = read_sigma(nlev, dl, du,Mobj)
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
if nargin >= 4
    conf.nlev = nlev;
    Mobj.sigma = dist;
    Mobj.siglev = zeros(Mobj.nVerts,conf.nlev);
    Mobj.siglevc = zeros(Mobj.nElems,conf.nlev);
    Mobj.siglayc = zeros(Mobj.nElems,conf.nlev-1);
    Mobj.siglev = repmat(Mobj.sigma,Mobj.nVerts,1);
    Mobj.siglay = Mobj.siglev(:,1:end-1) + (diff(Mobj.siglev,1,2)/2);
    for zz = 1:conf.nlev-1
        Mobj.siglevc(:, zz) = nodes2elems(Mobj.siglev(:, zz), Mobj);
        Mobj.siglayc(:, zz) = nodes2elems(Mobj.siglay(:, zz), Mobj);
    end
    % An extra level compared with layers.
    Mobj.siglevc(:, zz + 1) = nodes2elems(Mobj.siglev(:, zz + 1), Mobj);
    
    % Finally, make some [depth, sigma] arrays.
    Mobj.siglevz = repmat(Mobj.h, 1, conf.nlev) .* Mobj.siglev;
    Mobj.siglayz = repmat(Mobj.h, 1, conf.nlev-1) .* Mobj.siglay;
    if isfield(Mobj, 'hc')
        Mobj.siglevzc = repmat(Mobj.hc, 1, conf.nlev) .* Mobj.siglevc;
        Mobj.siglayzc= repmat(Mobj.hc, 1, conf.nlev-1) .* Mobj.siglayc;
    end
else
    Mobj = dist;
end
% generate sigma file 
% Save to the given file name.
if nargin==5
fout = fopen(sigma_file, 'wt');
assert(fout >= 0, 'Error opening sigma file: %s', sigma_file)
fprintf(fout, 'NUMBER OF SIGMA LEVELS = %d\n', nlev);
fprintf(fout, 'SIGMA COORDINATE TYPE = TANH\n');
fprintf(fout, 'DU = %4.1f\n', du);
fprintf(fout, 'DL = %4.1f\n', dl);
fprintf(fout,'\n');
fclose(fout);
end
return
