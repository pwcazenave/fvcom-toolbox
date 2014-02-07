function [Mobj] = read_gmsh_mesh(varargin)
% Read gmsh mesh files (version 2.2) into MATLAB mesh object.
%
% Mobj = function read_fvcom_mesh(varargin)
%
% DESCRIPTION:
%    Read gmsh msh file and bathymetry file. Store in a matlab mesh object.
%
% INPUT [keyword pairs]:
%   'msh'                   = gmsh msh file [e.g. casename.msh]
%   [optional] 'bath'       = gmsh bathymetry file [e.g. tst_dep.dat]
%   [optional] 'coordinate' = coordinate system [spherical; cartesian (default)]
%   [optional] 'addCoriolis' = calculate Coriolis param (f), requires [lon, lat]
%
% OUTPUT:
%    Mobj = MATLAB structure containing mesh data
%
% EXAMPLE USAGE
%    Mobj = read_gmsh_mesh('msh', 'casename.msh', 'bath', 'bathy.dat')
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%
%   2014-02-07 First version.
%
%==============================================================================

subname = 'read_gmsh_mesh';
global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

have_bath = false;
have_strings = false;

%--------------------------------------------------------------------------
% Create a blank mesh object
%--------------------------------------------------------------------------
Mobj = make_blank_mesh;
coordinate = 'cartesian';

%--------------------------------------------------------------------------
% Parse input arguments
%--------------------------------------------------------------------------

assert(mod(length(varargin), 2) == 0, 'incorrect usage of read_gmsh_mesh, use keyword pairs')

% Assume we have nothing sensible.
have_msh = false;
have_bath = false;
have_lonlat = false;
have_xy = false;

for i = 1:2:length(varargin) - 1
    keyword = lower(varargin{i});

    assert(ischar(keyword), 'incorrect usage of read_gmsh_mesh')

    switch keyword
        case 'msh'
            gmsh_msh = varargin{i + 1};
            have_msh = true;
        case 'bath'
            gmsh_bath = varargin{i + 1};
            have_bath = true;
        case 'coordinate'
            coord = varargin{i + 1};
            if strcmpi(coord, 'spherical')
                coordinate = 'spherical';
                have_lonlat = true;
            elseif strcmpi(coord, 'cartesian')
                coordinate = 'cartesian';
                have_xy = true;
            else
                warning('Unrecognised coordinate system (%s). Valid values are ''spherical'' and ''cartesian''.', coordinate)
            end
        case 'addcoriolis'
            val = varargin{i + 1};
            if val
                addCoriolis = true;
            else
                addCoriolis = false;
            end
        otherwise
            disp(varargin{i + 1})
            error('Can''t understand property: %s', varargin{i + 1});

    end
end

%--------------------------------------------------------------------------
% Read the mesh from the msh file
%--------------------------------------------------------------------------

fid = fopen(gmsh_msh, 'rt');
assert(fid >= 0, sprintf('file: %s does not exist\n',  gmsh_msh));

% Count number of elements and vertices
if ftbverbose
    fprintf('Reading from: %s\n', gmsh_msh)
end

% Read mesh type, taken from:
%
% http://www.geuz.org/pipermail/gmsh/attachments/20071002/642cb6c3/attachment.m
%

lin = fgetl(fid);
if strcmp(lin, '$NOD')
    fileformat = 1;
elseif strcmp(lin, '$MeshFormat')
    fileformat = 2;
    lin = fgetl(fid);
    if feof(fid)
        fprintf(sprintf('Syntax error (no version) in: %s\n', filename));
        fileformat = 0;
    end
    form = sscanf(lin, '%f %d %d');
    if form(1) < 2
        fprintf(sprintf('Unknown mesh format: %s\n', lin));
        fileformat = 0;
    end
    if form(2) ~= 0
        fprintf('Sorry, this program can only read ASCII format\n');
        fileformat = 0;
    end
    fgetl(fid); % this should be $EndMeshFormat
    if feof(fid)
        fprintf('Syntax error (no $EndMeshFormat) in: %s\n', filename);
        fileformat = 0;
    end
    lin = fgetl(fid); % this should be $Nodes
    if feof(fid)
        fprintf('Syntax error (no $Nodes) in: %s\n', filename);
        fileformat = 0;
    end
end

assert(logical(fileformat), 'Unrecognised mesh.')
% Read in the number of nodes
nVerts = str2double(fgetl(fid));

% Read the node positions.
C = textscan(fid, '%d %f %f %f', nVerts);
nid = C{1};
x = C{2};
y = C{3};
z = C{4};

if have_lonlat
    lon = x;
    lat = y;
else
    lon = zeros(size(x));
    lat = zeros(size(y));
end

% Now we should be at the end of the nodes and about to read the elements.
lin = fgetl(fid); lin = fgetl(fid); % need to skip one for some reason.
assert(strcmp(lin, '$EndNodes'), 'Improperly formatted msh file.')
lin = fgetl(fid);
assert(strcmp(lin, '$Elements'), 'Improperly formatted msh file.')
nElems = str2double(fgetl(fid));

% Preallocate some arrays of use. We'll trim the NaNs later.
tri12 = nan(nElems, 3);

% Read the element triangulation table. Format is:
%   ID, dim, n1, n2, n3[, n4]
% If dim is 1, then n4 is omitted, otherwise if dim is 2, then we have an
% extra link.
c1 = 1;
c2 = 1;
for i = 1:nElems
    lin = fgetl(fid);
    C = regexpi(lin, ' ', 'split');
    if str2double(C(2)) == 1
        tri12(c1, :) = str2double(C(end-2:end));
        c1 = c1 + 1;
    elseif str2double(C(2)) == 2
        tri12(c2, :) = str2double(C(end-2:end));
        c2 = c2 + 1;
%         else
%             warning('Unsupported number of dimensions (up to 2).')
    end
end
tri1 = tri12(1:c1, :);
tri2 = tri12(c1 + 1:c1 + c2, :);

% Trim the excess NaNs
tri2 = tri2(~isnan(tri2(:, 1)), :);

% Reset the number of elements to the actual number of 2D ones.
nElems = size(tri2, 1);

have_lonlat = false;
have_xy     = false;
if strcmpi(coordinate, 'spherical')
    lon = x;
    lat = y;
    % Why reset everything to zero here?
    %x = x * 0.0;
    %y = y * 0.0;
    have_lonlat = true;
    % Just do a double check on the coordinates to make sure we don't
    % actually have cartesian
    if max(lon) > 360
        warning('You''ve specified spherical coordinates, but your upper longitude value exceeds 360 degrees. Are you sure you have spherical data?')
    end
elseif strcmpi(coordinate, 'cartesian')
    have_xy = true;
else
    warning('Unrecognised coordinate system (%s). Valid values are ''spherical'' and ''cartesian''.', coordinate)
end

fclose(fid);

%--------------------------------------------------------------------------
% Read the topography from the bathymetry file
%--------------------------------------------------------------------------

if have_bath
    bath_range = max(h) - min(h);
    if have_bath || bath_range == 0
        fid = fopen(gmsh_bath, 'rt');
        if fid < 0
            error('File: %s does not exist', gmsh_bath);
        else
            if ftbverbose; fprintf('reading gmsh bathymetry from: %s\n', gmsh_bath); end
        end
        lin = fgetl(fid);
        lin = fgetl(fid);
        lin = fgetl(fid);
        C = textscan(fid, '%s %d', 1);
        nVerts_tmp = C{2};
        C = textscan(fid, '%s %d', 1);
        nElems_tmp = C{2};
        if (nVerts - nVerts_tmp) * (nElems - nElems_tmp) ~= 0
            fprintf('Dimensions of bathymetry file do not match msh file\n')
            fprintf('Bathymetry nVerts: %d\n', nVerts_tmp)
            fprintf('Bathymetry nElems: %d\n', nElems_tmp)
            error('Stopping...')
        end
        lin = fgetl(fid);
        lin = fgetl(fid);
        lin = fgetl(fid);
        lin = fgetl(fid); % extra one for the new format from SMS 10.1, I think
        C2 = textscan(fid, '%f', nVerts);
        h = C2{1};
        have_bath = true;

        clear C2
    elseif bath_range ~= 0
        have_bath = true;
    end

    % Make sure we have positive depths
    if sum(h > 0) < sum(h < 0)
        h = -h;
    end
end

%--------------------------------------------------------------------------
% Transfer to Mesh structure
%--------------------------------------------------------------------------

Mobj.nVerts  = nVerts;
Mobj.nElems  = nElems;
Mobj.nativeCoords = coordinate;

if have_lonlat
    Mobj.have_lonlat    = have_lonlat;
end
if have_xy
    Mobj.have_xy        = have_xy;
end
if have_bath
    Mobj.have_bath      = have_bath;
    Mobj.h              = h;
end
if have_strings
    Mobj.have_strings   = have_strings;
    Mobj.read_obc_nodes = read_obc_nodes;
end
if exist('addCoriolis', 'var') && addCoriolis
    Mobj.have_cor       = true;
end

Mobj.x            = x;
Mobj.y            = y;
Mobj.lon          = lon;
Mobj.lat          = lat;
Mobj.tri          = tri2;

if ftbverbose
  fprintf('end   : %s\n', subname)
end
