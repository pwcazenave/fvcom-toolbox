function [Mobj] = read_grid_mesh(varargin) 

% Read .grid mesh files into Matlab mesh object  
%
% [Mobj] = function read_grid_mesh(varargin)
%
% DESCRIPTION:
%    Read NOCL .grid file 
%    Store in a matlab mesh object 
%
% INPUT [keyword pairs]:  
%   'grid'                  = NOCL .grid file (e.g. UK_grid.grid)
%   [optional] 'coordinate' = coordinate system for output data [spherical; cartesian (default)]
%   [optional] 'input_coord' = coordinate system for input data [spherical; cartesian (default)]
%   [optional] 'project'    = generate (x,y) coordinates if input is (lon,lat) 
%                             generate (lon,lat) coordinates if input is (x,y)
%                            ['true' ; 'false' (default)], see myproject.m
%   [optional] 'zone'    = specify UTM zone for projection
%   [optional] 'addCoriolis' = calculate Coriolis param (f), requires [lon,lat]
%
% OUTPUT:
%    Mobj = matlab structure containing mesh data
%
% EXAMPLE USAGE
%    Mobj = read_grid_mesh('grid','UK_grid.grid','coordinate','spherical','addCoriolis','true')
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Karen Thurston (National Oceanography Centre Liverpool)
%
% Revision history (KJT)
%   2012-11-13 Adapted 'read_sms_mesh.m' to read NOCL .grid files
%   2012-11-19 Added input/output coordinate functionality
%   
%==============================================================================

subname = 'read_grid_mesh';
global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

have_grid = false;
have_bath = false;
have_lonlat = false;
have_xy = false;
userproject = false;
haveUTM = false;
addCoriolis = false;

%------------------------------------------------------------------------------
% Create a blank mesh object
%------------------------------------------------------------------------------
Mobj = make_blank_mesh();

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------

if(mod(length(varargin),2) ~= 0)
	error(['incorrect usage of ',subname,', use keyword pairs'])
end;


for i=1:2:length(varargin)-1
	keyword  = lower(varargin{i});
	if( ~ischar(keyword) )
		error(['incorrect usage of ',subname,', check keywords'])
	end;
	
	switch(keyword(1:3))
	
    case 'gri'
		grid_fname = varargin{i+1};
		have_grid = true;
	case 'coo'
		coord = varargin{i+1};
		if(coord(1:3)=='sph')
			coordinate = 'spherical';
		else
			coordinate = 'cartesian';
		end;
    case 'in_'
		in_coord = varargin{i+1};
		if(in_coord(1:3)=='sph')
			in_coordinate = 'spherical';
            have_lonlat = true;
		else
			in_coordinate = 'cartesian';
            have_xy = true;
		end;
	case 'pro'
		val = varargin{i+1};
		if( val )
			userproject = true;
		else
			userproject = false;
		end;
    case 'zon'
        fullzone = varargin{i+1};
        UTMzone = regexpi(fullzone,'\ ','split');
        UTMzone=str2double(char(UTMzone{1}(1)));
        haveUTM = true;
    case 'add'
		val = varargin{i+1};
		if( val )
			addCoriolis = true;
        else
			addCoriolis = false;
		end;
	otherwise
		error(['Can''t understand property:' char(varargin{i+1})]);
	end; %switch(keyword)
	
end;
		
%------------------------------------------------------------------------------
% Read the mesh from the .grid file
%------------------------------------------------------------------------------

fid = fopen(grid_fname,'r');
if(fid  < 0)
	error(['file: ' grid_fname ' does not exist']);
end;

% Count number of elements and vertices
if(ftbverbose);
  fprintf(['reading from: ' grid_fname '\n'])
end;
lin = fgetl(fid);   % ignore first line, it's the mesh name
c=textscan(fid,'%u %u',1);  % get nElems and nVerts
nElems = c{1};
nVerts = c{2};
clear c

if(ftbverbose); 
  fprintf('nVerts: %d\n',nVerts); 
  fprintf('nElems: %d\n',nElems); 
  fprintf('reading in connectivity and grid points\n')
end;

% allocate memory to hold mesh and connectivity
tri = zeros(nElems,3);
x   = zeros(nVerts,1);
y   = zeros(nVerts,1);
h   = zeros(nVerts,1);
lon = zeros(nVerts,1);
lat = zeros(nVerts,1);
ts  = zeros(nVerts,1);

c = textscan(fid, '%u %f %f %f ', nVerts);
x = c{2};
y = c{3};
h = c{4};
clear c

c = textscan(fid, '%u %u %u %u %u', nElems);
tri(:,1) = c{3};
tri(:,2) = c{4};
tri(:,3) = c{5};
clear c

% Make sure we have bathymetry
if sum(h)==0
    have_bath=false;
else
    have_bath=true;
end

% Make sure we have positive depths
if sum(h>0) < sum(h<0)
    h = -h;
end

% Build array of all the nodes in the open boundaries
c = textscan(fid, '%u %*[^\n]',1);
nOpenSeg = c{1};    % number of open boundary segments
clear c

lin=fgetl(fid); % skip the next line

c = textscan(fid, '%u %*[^\n]',1);
nOpenNodes = c{1};    % number of open boundary nodes
clear c

% Initialise SegCount variable
SegCount = [0,0];

for i=1:nOpenSeg    % for each open boundary segment
    c = textscan(fid, '%u %*[^\n]',1); % how many nodes in this segment?
    SegCount(1) = 1+SegCount(2);
    SegCount(2) = SegCount(1)+c{1}-1;
    clear c
    c = textscan(fid,'%u %*[^\n]',(SegCount(2)-SegCount(1)+1));    % get all the nodes in this segment
    allNodes{i} = c{1};
    clear c
end

if have_lonlat == true
	lon = x;
	lat = y;
	x = x*0.0;
	y = y*0.0;
	% Just do a double check on the coordinates to make sure we don't
    % actually have cartesian
    if max(lon)>360
        warning('You''ve specified spherical coordinates, but your upper longitude value exceeds 360 degrees. Are you sure you have spherical data?')
    end
else
	have_xy = true;
end;

%------------------------------------------------------------------------------
% Project if desired by user
%------------------------------------------------------------------------------

if(userproject)
    if (in_coordinate(1:3)=='car')
        fprintf('inverse projecting to get (lon,lat)\n')
        utmZones=cellfun(@(x) repmat(x,length(x),1),fullzone,'uni',false);
        [lon,lat] = utm2deg(x,y,utmZones{1});
        Mobj.have_lonlat = true;
    elseif (in_coordinate(1:3)=='sph')
        fprintf('forward projecting to get (x,y)\n')
        [x,y] = wgs2utm(lat,lon,UTMzone,'N');
        have_xy = true;
    end
end

%------------------------------------------------------------------------------
% Transfer to Mesh structure
%------------------------------------------------------------------------------

Mobj.nVerts  = nVerts;
Mobj.nElems  = nElems;
Mobj.nativeCoords = coordinate;

if have_lonlat==true
	Mobj.have_lonlat    = have_lonlat;
end;
if have_xy==true
	Mobj.have_xy        = have_xy;
end;

Mobj.have_bath      = have_bath;

Mobj.read_obc_nodes = allNodes;

Mobj.x            = x;
Mobj.y            = y;
Mobj.ts           = ts;
Mobj.lon          = lon;
Mobj.lat          = lat;
Mobj.h            = h;
Mobj.tri          = tri;

if addCoriolis==true
    Mobj.have_cor = true;
    Mobj = add_coriolis(Mobj,'uselatitude');
end

if(ftbverbose);
  fprintf(['end   : ' subname '\n'])
end;
fclose(fid);


