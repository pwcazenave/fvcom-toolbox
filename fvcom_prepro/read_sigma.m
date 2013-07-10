function Mobj = read_sigma(Mobj, sigmafile)
% Read an FVCOM sigma layers file and output z values into Mobj.
% 
% Mobj = read_sigma(Mobj, sigmafile)
% 
% DESCRIPTION:
%   Read a sigma file and calculate the sigma layer depths
% 
% INPUT:
%   Mobj:       Mesh object which must contain Mobj.h (depths).
%   sigmafile : Full path to an FVCOM sigma.dat file.
% 
% OUTPUT:
%   Mobj:       Mesh object with four new fields:
%               - siglayz and siglevz: contain depths of the sigma layers
%               and levels at each grid node.
%               - siglay and siglev: the sigma layer and levels in the
%               range 0 to -1.
% 
% EXAMPLE USAGE:
%   Mobj = read_sigma(Mobj, 'sigma.dat')
% 
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
% 
% Revision history
%   2013-01-08 Based on the code in show_sigma.m but instead of calculating
%   sigma layers along a user-defined line, the depths are calculated for
%   each node in the unstructured grid.
%   2013-01-10 Added two new outputs to the returned Mobj (siglay and
%   siglev). They're useful in write_FVCOM_tsobc.m.

subname = 'read_sigma';

global ftbverbose;
if ftbverbose
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end

fid = fopen(sigmafile,'r');
if(fid  < 0)
    error(['file: ' sigmafile ' does not exist']);
end

while ~feof(fid)
    line = fgetl(fid);
    if isempty(line) || strncmp(line, '!', 1) || ~ischar(line)
        continue
    end
    key = lower(line(1:3));
    C = strtrim(regexpi(line, '=', 'split'));
    switch key
        case 'num'
            nlev = str2double(C{2});
        case 'sig'
            sigtype = C{2};
        case 'du '
            du = str2double(C{2});
        case 'dl '
            dl = str2double(C{2});
        case 'min'
            min_constant_depth = str2double(C{2});
        case 'ku '
            ku = str2double(C{2});
        case 'kl '
            kl = str2double(C{2});
        case 'zku'
            s = str2double(regexp(C{2}, ' ', 'split'));
            zku = zeros(ku, 1);
            for i = 1:ku
                zku(i) = s(i);
            end
        case 'zkl'
            s = str2double(regexp(C{2}, ' ', 'split'));
            zkl = zeros(kl, 1);
            for i = 1:kl
                zkl(i) = s(i);
            end
    end
end

if ftbverbose
    fprintf('nlev %d\n',nlev)
    fprintf('sigtype %s\n',sigtype)
    fprintf('du %d\n',du)
    fprintf('dl %d\n',dl)
    fprintf('min_constant_depth %f\n',min_constant_depth)
    fprintf('ku %d\n',ku)
    fprintf('kl %d\n',kl)
    fprintf('zku %d\n',zku)
    fprintf('zkl %d\n',zkl)
end

% calculate the sigma distributions at each grid node
switch lower(sigtype)
    case 'generalized'
        z = sigma_gen(nlev, dl, du, kl, ku, zkl, zku, ...
            Mobj.z(i), min_constant_depth);
    case 'uniform'
        z = 0:-1/double(nlev-1):-1;
    otherwise
        error('Can''t do that sigtype')
end

% Create a siglay variable (i.e. midpoint in the sigma levels).
zlay = z(1:end-1) + (diff(z)/2);

Mobj.siglevz = repmat(Mobj.h, 1, nlev) .* repmat(z, Mobj.nVerts, 1);
Mobj.siglayz = repmat(Mobj.h, 1, nlev-1) .* repmat(zlay, Mobj.nVerts, 1);

% Add the sigma levels and layers to the Mobj.
Mobj.siglev = z;
Mobj.siglay = zlay;

if ftbverbose;
    fprintf(['end   : ' subname '\n'])
end