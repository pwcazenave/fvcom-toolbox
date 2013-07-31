function [Mobj] = read_sms_mesh(varargin) 
% Read sms mesh files into Matlab mesh object.
%
% [Mobj] = function read_fvcom_mesh(varargin)
%
% DESCRIPTION:
%    Read SMS 2dm file and bathymetry file 
%    Store in a matlab mesh object 
%
% INPUT [keyword pairs]:  
%   '2dm'                   = sms 2dm file [e.g. tst_grd.dat] 
%   [optional] 'bath'       = sms bathymetry file [e.g. tst_dep.dat] 
%   [optional] 'coordinate' = coordinate system [spherical; cartesian (default)]
%   [optional] 'project'    = generate (x,y) coordinates if input is (lon,lat) 
%                             generate (lon,lat) coordinates if input is (x,y)
%                            ['true' ; 'false' (default)], see myproject.m
%   [optional] 'addCoriolis' = calculate Coriolis param (f), requires [lon,lat]
%
% OUTPUT:
%    Mobj = matlab structure containing mesh data
%
% EXAMPLE USAGE
%    Mobj = read_sms_mesh('2dm','skagit.2dm','bath','bathy.dat')
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%
%   2012-06-20 Add support for reading nodestrings from SMS meshes.
%   2012-06-26 Added more resilient support for reading in SMS files.
%   2012-06-29 Further improved ability to read files with variable length
%   headers.
%   2013-07-31 Added some performance improvements to speed up loading mesh
%   files (from ~70s to ~30s on a 250,000 node grid). There's probably more
%   gains to be had by saving the values of tri, x, y and h when first
%   parsing the input file (lines 132-152). My brief testing would suggest
%   the overhead of converting from strings to doubles shouldn't be
%   underestimated.
%
%==============================================================================

subname = 'read_sms_mesh';
global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

userproject = false;
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

if mod(length(varargin), 2) ~= 0
	error('incorrect usage of read_sms_mesh, use keyword pairs')
end

for i = 1:2:length(varargin) - 1
	keyword = lower(varargin{i});

	if ~ischar(keyword)
		error('incorrect usage of read_sms_mesh')
    end

	switch keyword
        case '2dm'
            sms_2dm = varargin{i + 1};
            have_2dm = true;
        case 'bath'
            sms_bath = varargin{i + 1};
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
        case 'project'
            val = varargin{i + 1};
            if val
                userproject = true;
            else
                userproject = false;
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
% Read the mesh from the 2dm file
%--------------------------------------------------------------------------

fid = fopen(sms_2dm, 'rt');
if fid < 0
    error(['file: ' sms_2dm ' does not exist']);
end

% Count number of elements and vertices
if ftbverbose
    fprintf('reading from: %s\n', sms_2dm)
    fprintf('first pass to count number of nodes and vertices\n')
end

nElems = 0;
nVerts = 0;
nStrings = 0;
nHeader = 0;
StillReading = true;
while StillReading
	lin = fgetl(fid);
    if lin ~= -1 % EOF is -1
        switch lin(1:2)
            case 'E3'
                nElems = nElems + 1;
            case 'ND'
                nVerts = nVerts + 1;
            case 'NS'
                nStrings = nStrings + 1;
            case 'ME'
                nHeader = nHeader + 1;
            otherwise
                StillReading = false;
        end
    else
        % Got to EOF
        StillReading = false;
    end
end
fclose(fid);

fid = fopen(sms_2dm, 'rt');

if ftbverbose
  fprintf('nVerts: %d\n', nVerts); 
  fprintf('nElems: %d\n', nElems); 
  fprintf('reading in connectivity and grid points\n')
end

% Allocate memory to hold mesh and connectivity
tri = zeros(nElems,3);
lon = zeros(nVerts,1);
lat = zeros(nVerts,1);
ts  = zeros(nVerts,1);

% Skip the header
C = textscan(fid, '%s', nHeader + 1);

% Read the triangulation table
C = textscan(fid, '%s %d %d %d %d %d', nElems);
tri(:, 1) = C{3};
tri(:, 2) = C{4};
tri(:, 3) = C{5};

% Read in the nodes and interpolated depths
C = textscan(fid, '%s %d %f %f %f ', nVerts);
x = C{3};
y = C{4};
h = C{5};

% Check we don't have any NaNs anywhere
if max(isnan(x)) == 1
    error('%d NaNs in the x data', sum(isnan(x)))
end
if max(isnan(y)) == 1
    error('%d NaNs in the y data', sum(isnan(x)))
end
if max(isnan(h)) == 1
    error('%d NaNs in the h data', sum(isnan(x)))
end
if max(isnan(tri(:))) == 1
    error('%d NaNs in the h data', sum(isnan(tri(:))))
end

% Build array of all the nodes in the nodestrings.
C = textscan(fid, '%s %d %d %d %d %d %d %d %d %d %d', nStrings);
allNodes = cell2mat(C(2:end))';
allNodes(allNodes == 0) = [];

% Add a new field to Mobj with all the nodestrings as a cell array.
nodeStrings = find(allNodes < 0);
read_obc_nodes = cell(1, length(nodeStrings));
for nString = 1:sum(allNodes < 0)
    if nString == 1
        read_obc_nodes{nString} = abs(allNodes(1:nodeStrings(nString)));
    else
        read_obc_nodes{nString} = ...
            abs(allNodes(nodeStrings(nString - 1) + 1: ...
            nodeStrings(nString)));
    end
end

if nStrings > 0
    have_strings = true;
end

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

bath_range = max(h) - min(h);
if have_bath || bath_range == 0
	fid = fopen(sms_bath, 'rt');
	if fid < 0
		error('file: %s does not exist', sms_bath);
	else
		if ftbverbose; fprintf('reading sms bathymetry from: %s\n', sms_bath); end
    end
	lin = fgetl(fid);
	lin = fgetl(fid);
	lin = fgetl(fid);
	C = textscan(fid, '%s %d', 1);
	nVerts_tmp = C{2};
	C = textscan(fid, '%s %d', 1);
	nElems_tmp = C{2};
	if (nVerts - nVerts_tmp) * (nElems - nElems_tmp) ~= 0
		fprintf('dimensions of bathymetry file do not match 2dm file\n')
		fprintf('bathymetry nVerts: %d\n',nVerts_tmp)
		fprintf('bathymetry nElems: %d\n',nElems_tmp)
		error('stopping...')
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

%--------------------------------------------------------------------------
% Project if desired by user
%--------------------------------------------------------------------------
if userproject
    if strcmpi(coordinate, 'cartesian')
		fprintf('inverse projecting to get (lon,lat)\n')
		[lon, lat] = my_project(x, y, 'inverse');
		have_lonlat = true;
    elseif strcmpi(coordinate, 'spherical')
		fprintf('forward projecting to get (x,y)\n')
		[x, y] = my_project(lon, lat, 'forward');
		have_xy = true;
    else
        warning('Unrecognised coordinate system (%s). Valid values are ''spherical'' and ''cartesian''.', coordinate)
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
Mobj.ts           = ts;
Mobj.lon          = lon;
Mobj.lat          = lat;
Mobj.h            = h;
Mobj.tri          = tri;

if ftbverbose
  fprintf('end   : %s\n', subname)
end

fclose(fid);
