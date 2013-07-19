function [Mobj] = read_sms_mesh(varargin) 

% Read sms mesh files into Matlab mesh object  
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
%   
%==============================================================================

subname = 'read_sms_mesh';
global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

userproject = false;
have_bath = false;
have_strings = false;

%------------------------------------------------------------------------------
% Create a blank mesh object
%------------------------------------------------------------------------------
Mobj = make_blank_mesh();
coordinate = 'cartesian';

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------

if(mod(length(varargin),2) ~= 0)
	error('incorrect usage of read_sms_mesh, use keyword pairs')
end;


for i=1:2:length(varargin)-1
	keyword  = lower(varargin{i});
	if( ~ischar(keyword) )
		error('incorrect usage of read_sms_mesh')
	end;
	
	switch(keyword(1:3))
	
	case'2dm'
		sms_2dm = varargin{i+1};
		have_2dm = true;
	case 'bat'
		sms_bath = varargin{i+1};
		have_bath = true;
	case 'coo'
		coord = varargin{i+1};
		if(coord(1:3)=='sph')
			coordinate = 'spherical';
			have_lonlat = true;
		else
			coordinate = 'cartesian';
			have_xy    = true;
		end;
	case 'pro'
		val = varargin{i+1};
		if( val )
			userproject = true;
		else
			userproject = false;
		end;
	case 'add'
		val = varargin{i+1};
		if( val )
			addCoriolis = true;
        else
			addCoriolis = false;
		end;
	otherwise
		error(['Can''t understand property:' varargin{i+1}]);
	end; %switch(keyword)
	
end;
		
%------------------------------------------------------------------------------
% Read the mesh from the 2dm file
%------------------------------------------------------------------------------


fid = fopen(sms_2dm,'r');
if(fid  < 0)
	error(['file: ' sms_2dm ' does not exist']);
end;

% Count number of elements and vertices
if(ftbverbose);
  fprintf(['reading from: ' sms_2dm '\n'])
  fprintf('first pass to count number of nodes and vertices\n')
end;
nElems = 0;
nVerts = 0;
nStrings = 0;
nHeader = 0;
StillReading = true;
while StillReading
	lin = fgetl(fid);
    if lin ~= -1 % EOF is -1
        if(lin(1:2) == 'E3')
            nElems = nElems + 1;
        elseif(lin(1:2) == 'ND')
            nVerts = nVerts + 1;
        elseif(lin(1:2) == 'NS')
            nStrings = nStrings + 1;
        elseif(lin(1:2) == 'ME')
            % Header values
            nHeader = nHeader + 1;
        else
            StillReading = false;
        end
    else
        % Got to EOF
        StillReading = false;
    end
end
fclose(fid);
fid = fopen(sms_2dm,'r');

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

validObs = 1;

% skip the header
for i=1:nHeader
    if i==1
        C = textscan(fid,'%s',1);
    else
        C = textscan(fid,'%s',2);
    end
end
clear C

for i=1:nElems
    C = textscan(fid, '%s %d %d %d %d %d', 1);
    % Check we have valid data.
    if strcmp(C{1},'E3T')
        tri(validObs,1) = C{3};
        tri(validObs,2) = C{4};
        tri(validObs,3) = C{5};
        validObs = validObs + 1;
    end
end

for i=1:nVerts 
	C = textscan(fid, '%s %d %f %f %f ', 1);
	x(i) = C{3};
	y(i) = C{4};
    h(i) = C{5};
end

% Build array of all the nodes in the nodestrings.
for i=1:nStrings
    C = textscan(fid, '%s %d %d %d %d %d %d %d %d %d %d', 1);
    C2 = cell2mat(C(2:end));
    C2(C2==0) = [];
    if i == 1
        allNodes = C2;
    else
        allNodes = [allNodes, C2];
    end
end

% Add a new field to Mobj with all the nodestrings as a cell array.
nodeStrings = find(allNodes<0);
for nString=1:sum(allNodes<0)
    if nString == 1
        read_obc_nodes{nString} = abs(allNodes(1:nodeStrings(nString)));
    else
        read_obc_nodes{nString} = abs(allNodes(nodeStrings(nString-1)+1:nodeStrings(nString)));
    end
end
% Just add them all to the obc_nodes
% obc_nodes = double(abs(allNodes));

if nStrings > 0
    have_strings = true;
end

have_lonlat = false;
have_xy     = false;
if(coordinate(1:5) == 'spher')
	lon = x;
	lat = y;
	x = x*0.0;
	y = y*0.0;
	have_lonlat = true;
    % Just do a double check on the coordinates to make sure we don't
    % actually have cartesian
    if max(lon)>360
        warning('You''ve specified spherical coordinates, but your upper longitude value exceeds 360 degrees. Are you sure you have spherical data?')
    end
else
	have_xy = true;
end;


%------------------------------------------------------------------------------
% Read the topography from the bathymetry file
%------------------------------------------------------------------------------

bath_range = max(h)-min(h);
if(have_bath) || bath_range == 0
	fid = fopen(sms_bath,'r');
	if(fid  < 0)
		error(['file: ' sms_bath ' does not exist']);
	else
		if(ftbverbose); fprintf('reading sms bathymetry from: %s\n',sms_bath); end;
	end;
	lin=fgetl(fid);
	lin=fgetl(fid);
	lin=fgetl(fid);
	C = textscan(fid, '%s %d', 1);
	nVerts_tmp = C{2};
	C = textscan(fid, '%s %d', 1);
	nElems_tmp = C{2};
	if( (nVerts-nVerts_tmp)*(nElems-nElems_tmp) ~= 0)
		fprintf('dimensions of bathymetry file do not match 2dm file\n')
		fprintf('bathymetry nVerts: %d\n',nVerts_tmp)
		fprintf('bathymetry nElems: %d\n',nElems_tmp)
		error('stopping...')
	end;
	lin=fgetl(fid);
	lin=fgetl(fid);
	lin=fgetl(fid);
    lin=fgetl(fid); % extra one for the new format from SMS 10.1, I think
	for i=1:nVerts
	  C = textscan(fid, '%f', 1);
	  h(i) = C{1}; 
	end;
	have_bath = true;
elseif bath_range ~= 0
    have_bath = true;
end;
% Make sure we have positive depths
if sum(h>0) < sum(h<0)
    h = -h;
end

%------------------------------------------------------------------------------
% Project if desired by user
%------------------------------------------------------------------------------

if(userproject)
	if(coordinate(1:3) == 'car')
		fprintf('inverse projecting to get (lon,lat)\n')
		[lon,lat] = my_project(x,y,'inverse');
		have_lonlat = true;
	else
		fprintf('forward projecting to get (x,y)\n')
		[x,y] = my_project(lon,lat,'forward');
		have_xy = true;
	end;
end;

%------------------------------------------------------------------------------
% Transfer to Mesh structure
%------------------------------------------------------------------------------

Mobj.nVerts  = nVerts;
Mobj.nElems  = nElems;
Mobj.nativeCoords = coordinate;

if(have_lonlat)
	Mobj.have_lonlat    = have_lonlat;
end;
if(have_xy)
	Mobj.have_xy        = have_xy;
end;
if(have_bath)
	Mobj.have_bath      = have_bath;
end;
if(have_strings)
    Mobj.have_strings   = have_strings;
    Mobj.read_obc_nodes = read_obc_nodes;
end
if exist('addCoriolis','var') && addCoriolis
    Mobj.have_cor       = true;
end

Mobj.x            = x;
Mobj.y            = y;
Mobj.ts           = ts;
Mobj.lon          = lon;
Mobj.lat          = lat;
Mobj.h            = h;
Mobj.tri          = tri;

if(ftbverbose);
  fprintf(['end   : ' subname '\n'])
end;
fclose(fid);


