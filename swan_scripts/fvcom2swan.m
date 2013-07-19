function fvcom2swan(fvcom_mesh,fvcom_bathy,fvcom_obc,prefix,PlotMesh)
% Convert fvcom grid and bathymetry file to unstructured SWAN format 
%
% function fvcom2swan(fvcom_mesh,fvcom_bathy,fvcom_obc,prefix)
%
% DESCRIPTION:
%    convert FVCOM mesh and bathymetry to SWAN format mesh and bathymetry
%
% INPUT 
%   fvcom_mesh  = FVCOM 3.x grid file 
%   fvcom_bathy = FVCOM 3.x bathymetry file
%   fvcom_obc   = FVCOM 3.x open boundary file
%   prefix      = prefix for naming SWAN grid and depth files   
%   PlotMesh    = [true,false] plot the resulting mesh          
%
% OUTPUT:
%    prefix.bot  = swan bathymetry file
%    prefix.node = swan vertex file
%    prefix.ele  = swan connectivity file
%
% EXAMPLE USAGE
%    fvcom2swan('tst_grd.dat','tst_dep.dat','tst_obc.dat','tst',true) 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

subname = 'fvcom2swan';
fprintf('\n')
fprintf(['begin : ' subname '\n'])

% setup SWAN output files
swan_bathy  = [prefix,'.bot'];
swan_node   = [prefix,'.node'];
swan_ele    = [prefix,'.ele'];

%------------------------------------------------------------------------------
% read in the FVCOM bathymetry data
%------------------------------------------------------------------------------
fid = fopen(fvcom_bathy,'r');
if(fid  < 0)
  error(['file: ' fvcom_bathy ' does not exist']);
end;
C = textscan(fid, '%s %s %s %d', 1);
Nverts = C{4};
h = zeros(Nverts,1);
fprintf('reading bathymetry file\n');
fprintf('# nodes %d\n',Nverts);
for i=1:Nverts
  C = textscan(fid, '%f %f %f', 1);
  h(i) = C{3};
end;
fprintf('min depth %f max depth %f\n',min(h),max(h));
fprintf('bathymetry reading complete\n');
fclose(fid);

%------------------------------------------------------------------------------
% read in the FVCOM open boundary node data
%------------------------------------------------------------------------------
fid = fopen(fvcom_obc,'r');
if(fid  < 0)
  error(['file: ' fvcom_obc ' does not exist']);
end;
C = textscan(fid, '%s %s %s %s %d', 1);
Nobcs = C{5};
obc_nodes = zeros(Nobcs,1);
fprintf('reading obc file\n');
fprintf('# nodes %d\n',Nobcs);
for i=1:Nobcs
  C = textscan(fid, '%d %d %d', 1);
  obc_nodes(i) = C{2};
end;

fprintf('obc reading complete\n');
fclose(fid);

%----------------------------------------------------
% read in the fvcom connectivity and vertices 
%----------------------------------------------------
fid = fopen(fvcom_mesh,'r');
if(fid  < 0)
  error(['file: ' fvcom_mesh ' does not exist']);
end;
C = textscan(fid, '%s %s %s %d', 1); Nverts = C{4};
C = textscan(fid, '%s %s %s %d', 1); Nelems = C{4};
tri = zeros(Nelems,3); 
x   = zeros(Nverts,2);
fprintf('reading mesh file\n');
fprintf('# nodes %d\n',Nverts);
fprintf('# elems %d\n',Nelems);
for i=1:Nelems
  C = textscan(fid,' %d %d %d %d %d\n',1);
  tri(i,1) = C{2};  tri(i,2) = C{3}; tri(i,3) = C{4};
end;
for i=1:Nverts 
  C = textscan(fid, '%d %f %f %f', 1);
  x(i,1) = C{2};
  x(i,2) = C{3};
end;
fprintf('mesh read in\n');
fclose(fid);

%----------------------------------------------------
% mark nodes on the boundary 
%----------------------------------------------------
bnodes = zeros(Nverts,1);
cells = zeros(Nverts,10);
cellcnt = zeros(Nverts,1);
nbe = zeros(Nelems,3);

for i = 1:Nelems
	n1 = tri(i,1) ; cellcnt(n1) = cellcnt(n1) + 1;
	n2 = tri(i,2) ; cellcnt(n2) = cellcnt(n2) + 1;
	n3 = tri(i,3) ; cellcnt(n3) = cellcnt(n3) + 1;
	cells(tri(i,1),cellcnt(n1)) = i;
	cells(tri(i,2),cellcnt(n2)) = i;
	cells(tri(i,3),cellcnt(n3)) = i;
end;

if(max(cellcnt) > 10)
  error('increase cells array')
end;

for i = 1:Nelems
	n1 = tri(i,1); n2 = tri(i,2); n3 = tri(i,3);
	for j1 = 1:cellcnt(n1)
		for j2 = 1:cellcnt(n2)
			if((cells(n1,j1) == cells(n2,j2)) & cells(n1,j1) ~= i); nbe(i,3) = cells(n1,j1); end;
		end;
	end;
	for j2 = 1:cellcnt(n2)
		for j3 = 1:cellcnt(n3)
			if((cells(n2,j2) == cells(n3,j3)) & cells(n2,j2) ~= i); nbe(i,1) = cells(n2,j2); end;
		end;
	end;
	for j1 = 1:cellcnt(n1)
		for j3 = 1:cellcnt(n3)
 			if((cells(n1,j1) == cells(n3,j3)) & cells(n1,j1) ~= i); nbe(i,2) = cells(n3,j3); end;
		end;
	end;
end;


check = Nelems/10;
for i=1:Nelems
   n1 = tri(i,1); n2 = tri(i,2); n3 = tri(i,3);
   if(nbe(i,1) == 0)
     bnodes(n2) = 1; bnodes(n3) = 1;
   elseif(nbe(i,2) == 0)
     bnodes(n3) = 1; bnodes(n1) = 1; 
   elseif(nbe(i,3) == 0)
     bnodes(n1) = 1; bnodes(n2) = 1; 
   end;
   if(mod(i,check)==0); fprintf('bnodes: completed %f percent \n',100*i/Nelems); end;
end;

%----------------------------------------------------
% Mark open boundary with value of 2 
%----------------------------------------------------
if(Nobcs > 1)
   bnodes(obc_nodes) = 2;
end;
%----------------------------------------------------
% dump swan bathymetry file 
%----------------------------------------------------
fid = fopen(swan_bathy,'w');
for i=1:Nverts
	fprintf(fid,'%f\n',h(i));
end;
fclose(fid);


%----------------------------------------------------
% dump swan node file 
%----------------------------------------------------
fid = fopen(swan_node,'w');
fprintf(fid,'%d 2 0 1\n',Nverts);
for i=1:Nverts
	fprintf(fid,'%d %f %f %d\n',i,x(i,1),x(i,2),bnodes(i)); 
end;
fclose(fid);


%----------------------------------------------------
% dump swan connectivty file 
%----------------------------------------------------
fid = fopen(swan_ele,'w');
fprintf(fid,'%d 3 0\n',Nelems);
for i=1:Nelems 
	fprintf(fid,'%d %d %d %d\n',i,tri(i,1:3)); 
end;
fclose(fid);

%----------------------------------------------------
% plot mesh from swan files to check 
%----------------------------------------------------
if(PlotMesh)
	plot_swan_mesh(swan_bathy,swan_node,swan_ele)
end;

fprintf(['end   : ' subname '\n'])

