
% create a list of the nesting boundary nodes
% subname = 'create_nesting_nodes'
%
% DESCRIPTION:
%    READ in the following data:
%      	- FVCOM OBC file  
%    	- FVCOM Grid file (connectivity + nodes)
%    Process data:
%       - Extract all nodes belonging to elements at the boundary 
%    Write output
%       - Write Node nesting file 
%
% INPUT :  
%   FVCOM OBC file and Grid file
%
% OUTPUT:
%   Node nesting file
%
% EXAMPLE USAGE
%    create_nesting_nodes
%
% Author(s):  
%    Hakeem Johnson (CH2MHILL, Warrington, UK)
%    Darren Price(CH2MHILL, Warrington, UK)
%
% Revision history
%   Jan 2014: Beta version including the correct nesting format
%==============================================================================

clear;
subname = 'create_nesting_nodes';

global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end

%---------------------------------------------------------
% user specify input data & parameters ...
%---------------------------------------------------------
% work & data folders ...
workDir='\\hand-fs-01\maritime\Projects\Scottish Waters\Calcs\Models\WLLS\mesh10\nesting\';
dataDir='\\hand-fs-01\maritime\Projects\Scottish Waters\Calcs\Models\WLLS\mesh10\mesh\';
basename = 'WLLS_v3smth';
% input file names
meshFile=[basename '_grd.dat'];
obcFile=[basename '_obc.dat'];

% output file
NestNodesFile = [basename '_node_nest.dat'];

%---------------------------------------------------------
% set full file names ...
%---------------------------------------------------------
meshFile     = [dataDir meshFile];
obcFile      = [dataDir obcFile];
NestNodesFile     = [workDir NestNodesFile]; 
checkfile(meshFile);
checkfile(obcFile);

%---------------------------------------------------------
% Read data into matlab 
% 1) grid mesh; 2)boundary nodes; ...
%---------------------------------------------------------
% get fvcom grid file in as mesh object ...
FV_Mesh = read_fvcom_mesh_lonlat(meshFile);		% ... from matlab toolbox 

% get wave nesting nodes and associated lat/lon coord file
FV_OBC = get_HD_OBC(obcFile);

%---------------------------------------------------------------------
% get elements at boundary & element centres
%---------------------------------------------------------------------
% get elements at boundary
T = FV_Mesh.tri;
X1 = FV_Mesh.lon;
Y1 = FV_Mesh.lat;
P = [X1,Y1];
TR = triangulation(T,P);
vi = FV_OBC.nnodesID;
ti = vertexAttachments(TR,vi);
% arrange as column data; since elements overlap, get unique elements
temp1 = [ti{:}]';
% bdcell = unique(temp1,'stable');  
bdcell = unique(temp1);  
nCells = numel(bdcell);

nElems = nCells;
bdElem = bdcell;

% get element nodes & arrange in required format
nv= zeros(3*nElems,1);
k = 1;
for i=1:nElems
    ielem = bdElem(i);
    nv(k)   = FV_Mesh.tri(ielem,1);
    nv(k+1) = FV_Mesh.tri(ielem,2);
    nv(k+2) = FV_Mesh.tri(ielem,3);
    k = k+3;
end

% I don't think preserving the order is the problem - remove next line
% nestNodeID= unique(nv,'stable');  
nestNodeID= unique(nv);
nestNodes = numel(nestNodeID);
out = ones(nestNodes,3);
for i=1:nestNodes
    out(i,1) = i;
    out(i,2) = nestNodeID(i);
end

%------------------------------------------------------------------------------
% Dump the file
%------------------------------------------------------------------------------

filename = NestNodesFile;
if(ftbverbose); fprintf('writing FVCOM obc %s\n',filename); end;

%------------------------------------------------------------------------------
% Check output file
%------------------------------------------------------------------------------
fid = fopen(filename,'w');
fprintf(fid,'Nest_Node Number = %12d\n',nestNodes);
fprintf(fid,'No.           Node_ID            Node_Type\n');
for i=1:nestNodes
	fprintf(fid,'%12d %12d %12d\n',out(i,:));
end

fclose(fid);

disp('finished creating nesting file');