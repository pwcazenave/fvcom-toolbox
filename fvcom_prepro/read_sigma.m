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
%               - siglayzc and siglevzc: contain depths of the sigma layers
%               and levels at each element centre.
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
%   2013-04-23 Add support for geometric sigma distributions as well as
%   slightly more robust reading of and checks on the parameters in the
%   input file. Also changed the way the uniform distribution is calculated
%   (by using a P_SIGMA value of 1 and the sigma_geo.m function rather than
%   fiddling around with ranges, although the output is the same).
%   2014-04-28 Add the sigma levels for the element centres in addition to
%   the element nodes.

subname = 'read_sigma';

global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
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

    % Clean up the input string to make matching a bit easier (trim
    % whitespace and remove duplicate spaces in the keywords).
    C = strtrim(regexpi(regexprep(line, '\s+', ' '), '=', 'split'));

    switch lower(C{1})
        case 'number of sigma levels'
            nlev = str2double(C{2});
        case 'sigma coordinate type'
            sigtype = C{2};
        case 'sigma power'
            sigpow = str2double(C{2});
        case 'du'
            du = str2double(C{2});
        case 'dl'
            dl = str2double(C{2});
        case 'min constant depth'
            min_constant_depth = str2double(C{2});
        case 'ku'
            ku = str2double(C{2});
        case 'kl'
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

% Do some checks if we've got uniform or generalised coordinates to make
% sure the input is correct.
if strcmpi(sigtype, 'GENERALIZED') || strcmpi(sigtype, 'UNIFORM')
    if numel(zku) ~= ku
        warning('Number of zku values does not match the number specified in ku')
    end
    if numel(zkl) ~= kl
        warning('Number of zkl values does not match the number specified in kl')
    end
end

if ftbverbose
    % Should be present in all sigma files.
    fprintf('nlev\t%d\n', nlev)
    fprintf('sigtype\t%s\n', sigtype)

    % Only present in geometric sigma files.
    if strcmpi(sigtype, 'GEOMETRIC')
        fprintf('sigpow\t%d\n', sigpow)
    end

    % Only in the generalised or uniform sigma files.
    if strcmpi(sigtype, 'GENERALIZED') || strcmpi(sigtype, 'UNIFORM')
        fprintf('du\t%d\n', du)
        fprintf('dl\t%d\n', dl)
        fprintf('min_constant_depth\t%f\n', min_constant_depth)
        fprintf('ku\t%d\n', ku)
        fprintf('kl\t%d\n', kl)
        fprintf('zku\t%d\n', zku)
        fprintf('zkl\t%d\n', zkl)
    end
end

% Calculate the sigma distributions at each grid node
switch lower(sigtype)
    case 'generalized'
        z = sigma_gen(nlev, dl, du, kl, ku, zkl, zku, ...
            Mobj.z(i), min_constant_depth);
    case 'uniform'
        z = sigma_geo(nlev, 1);
        % z = 0:-1/double(nlev-1):-1;
    case 'geometric'
        z = sigma_geo(nlev, sigpow);
    otherwise
        error('Don''t know sigtype %s (is it supported?)', sigtype)
end

% Create a siglay variable (i.e. midpoint in the sigma levels).
zlay = z(1:end-1) + (diff(z)/2);

% Create a depth array for the element centres.
hc = nodes2elems(Mobj.h, Mobj);

Mobj.siglevz = repmat(Mobj.h, 1, nlev) .* repmat(z, Mobj.nVerts, 1);
Mobj.siglayz = repmat(Mobj.h, 1, nlev-1) .* repmat(zlay, Mobj.nVerts, 1);
Mobj.siglevzc = repmat(hc, 1, nlev) .* repmat(z, Mobj.nElems, 1);
Mobj.siglayzc = repmat(hc, 1, nlev-1) .* repmat(zlay, Mobj.nElems, 1);

% Add the sigma levels and layers to the Mobj.
Mobj.siglev = z;
Mobj.siglay = zlay;

if ftbverbose;
    fprintf('end   : %s\n', subname)
end
