function [h] = read_fvcom_bath(bathfile)

% Read fvcom bathymetry file
%
% [h] = function read_fvcom_bath(bathfile)
%
% DESCRIPTION:
%    Read FVCOM Bathymetry file
%
% INPUT:
%   bathfile  = fvcom bathymetry file
%
% OUTPUT:
%    h = bathymetry vector
%
% EXAMPLE USAGE
%    Mobj = read_fvcom_bath('tst_dep.dat')
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Rory O'Hara Murray (Marine Scotland Science)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2014-11-19 Remove loops to speed up reading in the file.
%    2017-04-11 Minor clean up of the code.
%
%==============================================================================

global ftbverbose
[~, subname] = fileparts(mfilename('fullpath'));
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

%------------------------------------------------------------------------------
% read in the FVCOM bathymetry data
%------------------------------------------------------------------------------
fid = fopen(bathfile, 'r');
assert(fid > 0, 'file: %s does not exist', bathfile);

C = textscan(fid, '%s %s %s %d', 1);
Nverts = C{4};
if ftbverbose
    fprintf('reading bathymetry file\n');
    fprintf('# nodes %d\n',Nverts);
end
C = textscan(fid,' %f %f %f',Nverts);
h = C{3};
fclose(fid);

if ftbverbose
    fprintf('min depth %f max depth %f\n',min(h),max(h));
    fprintf('bathymetry reading complete\n');
end

if ftbverbose
    fprintf('end   : %s\n', subname)
end
