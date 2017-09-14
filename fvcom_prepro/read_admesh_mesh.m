function [Mobj] = read_admesh_mesh(varargin)
% Read admesh mesh files (version 2.2) into MATLAB mesh object.
%
% Mobj = function read_fvcom_mesh(varargin)
%
% DESCRIPTION:
%    Read admesh 14 file and bathymetry file. Store in a matlab mesh object.
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
%    Mobj = read_admesh_mesh('msh', 'casename.14', 'coordinate', 'spherical')
%
% Author(s):
%    Ricardo Torres (Plymouth Marine Laboratory) based on read_gmsh_mesh
%
% Revision history
%
%   2016-06-22 First version.
%
%==============================================================================

subname = 'read_admesh_mesh';
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

assert(mod(length(varargin), 2) == 0, 'incorrect usage of read_admesh_mesh, use keyword pairs')

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

% Read mesh type, written from srcatch:
%
% http://www.geuz.org/pipermail/gmsh/attachments/20071002/642cb6c3/attachment.m
%

lin = fgetl(fid);
% first line is title of mesh
title_str = lin;
lin = fgetl(fid);
% next line is mesh dimensions
form = sscanf(lin, '%u %u');
nElems = form(1);
nVerts = form(2);

% 
% if strcmp(lin, '$NOD')
%     fileformat = 1;
% elseif strcmp(lin, '$MeshFormat')
%     fileformat = 2;
%     lin = fgetl(fid);
%     if feof(fid)
%         fprintf(sprintf('Syntax error (no version) in: %s\n', filename));
%         fileformat = 0;
%     end
%     form = sscanf(lin, '%f %d %d');
%     if form(1) < 2
%         fprintf(sprintf('Unknown mesh format: %s\n', lin));
%         fileformat = 0;
%     end
%     if form(2) ~= 0
%         fprintf('Sorry, this program can only read ASCII format\n');
%         fileformat = 0;
%     end
%     fgetl(fid); % this should be $EndMeshFormat
%     if feof(fid)
%         fprintf('Syntax error (no $EndMeshFormat) in: %s\n', filename);
%         fileformat = 0;
%     end
%     lin = fgetl(fid); % this should be $Nodes
%     if feof(fid)
%         fprintf('Syntax error (no $Nodes) in: %s\n', filename);
%         fileformat = 0;
%     end
% end
% 
% assert(logical(fileformat), 'Unrecognised mesh.')
% Read in the number of nodes
% nVerts = str2double(fgetl(fid));

% Read the node positions and depths.
C = textscan(fid, '%d %f %f %f', nVerts);
nid = C{1};
x = C{2};
y = C{3};
h = C{4};

if have_lonlat
    lon = x;
    lat = y;
else
    lon = zeros(size(x));
    lat = zeros(size(y));
end

% Now we should be at the end of the nodes and about to read the elements.

% Read the element triangulation table. Format is:
%   ID, dim, n1, n2, n3
% Read the Element positions.
C = textscan(fid, '%d %d %d %d %d', nElems);

tri1 = [C{3},C{4},C{5}];

have_lonlat = false;
have_xy     = false;
if strcmpi(coordinate, 'spherical')
    lon = x;
    lat = y;
    % Why reset everything to zero here? Because we don't have cartesian
    % ... what we read before were lat and lon
    x = x * 0.0;
    y = y * 0.0;
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


%--------------------------------------------------------------------------
% check the topography from the mesh file... No separate file for
% bathymetry here
%--------------------------------------------------------------------------

    bath_range = max(h) - min(h);
    if have_bath || bath_range == 0
    elseif bath_range ~= 0
        have_bath = true;
    end
    % Make sure we have positive depths
    if sum(h > 0) < sum(h < 0)
        h = -h;
    end
%--------------------------------------------------------------------------
% Read the nodestrings 
%--------------------------------------------------------------------------
lin = fgetl(fid); % need to skip one line for some reason...
lin = fgetl(fid);
if feof(fid)
    fprintf(sprintf('No Open Boundary strings found in: %s\n', filename));
    have_strings=0;
else
    ind0=strfind(lin, '=');
    str_case = strtrim(lin(ind0+1:end));
    if strcmp(str_case, 'Number of open boundaries');
        nBoundary = sscanf(lin(1:ind0-1), '%d');
        have_strings=1;
    else
        fprintf(sprintf('Unknown mesh file format: %s\n', lin));
    end
end

if have_strings
    % read total number of open boundary nodes
     lin = fgetl(fid);
    ind0=strfind(lin, '=');
    str_case = strtrim(lin(ind0+1:end));
    if strcmp(str_case, 'Total number of open boundary nodes');
        nObc = sscanf(lin(1:ind0-1), '%d');
    else
        fprintf(sprintf('Unknown mesh file format: %s\n', lin));
    end
    % read number of open boundary nodes for each boundary
    for nn=1:nBoundary
    lin = fgetl(fid);
    ind0=strfind(lin, '=');
    str_case = strtrim(lin(ind0+1:end));
    display([' Reading  ',str_case ]);
        tmpnObc = sscanf(lin(1:ind0-1), '%u');
        boundary_names{nn} = num2str(nn);
        read_obc_nodes(nn)= textscan(fid,'%d',tmpnObc);
        fprintf(sprintf('Read Nodes in boundary: %s\n', boundary_names{nn}));
        lin = fgetl(fid); % no idea why we need this additional one...

    end
end
fclose(fid);


%--------------------------------------------------------------------------
% Transfer to Mesh structure
%--------------------------------------------------------------------------

Mobj.nVerts  = nVerts;
Mobj.nElems  = nElems;
Mobj.nativeCoords = coordinate;

if have_lonlat
    Mobj.have_lonlat    = have_lonlat;
Mobj.lon          = lon;
Mobj.lat          = lat;
end
if have_xy
    Mobj.have_xy        = have_xy;
Mobj.x            = x;
Mobj.y            = y;
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

Mobj.tri          = tri1;

if ftbverbose
  fprintf('end   : %s\n', subname)
end
