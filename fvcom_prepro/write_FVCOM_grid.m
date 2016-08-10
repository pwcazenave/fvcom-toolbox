function write_FVCOM_grid(Mobj,filename)
% Write grid and connectivity to FVCOM format grid file
%
% function write_FVCOM_grid(Mobj,filename)
%
% DESCRIPTION:
%    Generate an ascii FVCOM 3.x format gridfile from Mesh object
%
% INPUT
%   Mobj     = Mesh object with fields:
%       - nativeCoords - string of the native coordinates (cartesian or
%       spherical).
%       - x, y, lon, lat - node coordinates for either cartesian or
%       spherical.
%       - nVerts - number of nodes.
%       - nElems - number of elements.
%       - tri - grid connectivity table.
%       - h - water depth at the nodes.
%   filename = FVCOM grid file name
%
% OUTPUT:
%    FVCOM grid file: filename
%
% EXAMPLE USAGE
%    write_FVCOM_grid(Mobj, 'tst_grd.dat')
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Rory O'Hara Murray (Marine Scotland Science)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2014-10-07 Removed loops to speed up writing the file
%    2016-06-14 Minor code tweaks to fix MATLAB warnings. Add check that we
%    successfully opened the output file.
%==========================================================================

[~, subname] = fileparts(mfilename('fullpath'));

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname);
end

%--------------------------------------------------------------------------
% Parse input arguments
%--------------------------------------------------------------------------
if ~exist('Mobj', 'var') || ~exist('filename', 'var')
    error('arguments to write_FVCOM_grid are incorrect')
end

%--------------------------------------------------------------------------
% Dump the file
%------------------------------------------------------------------------------
if strcmpi(Mobj.nativeCoords, 'cartesian')
    x = Mobj.x;
    y = Mobj.y;
else
    x = Mobj.lon;
    y = Mobj.lat;
end
if ftbverbose
    fprintf('writing FVCOM gridfile %s\n', filename)
end
fid = fopen(filename,'w');
assert(fid > 0, 'Error opening output file %s', filename)
fprintf(fid, 'Node Number = %d\n', Mobj.nVerts);
fprintf(fid, 'Cell Number = %d\n', Mobj.nElems);
fprintf(fid, '%d %d %d %d %d\n', [(1:Mobj.nElems)', Mobj.tri, (1:Mobj.nElems)']');
fprintf(fid, '%d %f %f %f\n', [(1:Mobj.nVerts)', x, y, Mobj.h]');
fclose(fid);

if ftbverbose
    fprintf('end   : %s\n', subname)
end

