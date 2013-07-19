function write_dtank(fname,Mobj,geographic)
% Dump mesh to datatank file 
%
% function dump_dtascii(fname,Mobj) 
%
% DESCRIPTION:
%    Dump mesh, open boundary nodes, river nodes, and probe nodes to DTascii file 
%
% INPUT 
%   fname = datatank file name
%   Mobj  = mesh object
%
% OUTPUT:
%    fname:  ascii datatank file 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history

% parse input args
cartesian = true;
if(exist('geographic'))
  cartesian = false;
end;

if(~cartesian)
  % open datatank file
  fid = fopen(fname,'w');
  fprintf(fid,'DataTank Binary File\n');
  fprintf(fid,'Bathymetry_Geog\n');
  fprintf(fid,'TriangularMesh2D\n');
  
  %Grid
  fprintf(fid,'Points\n');
  fprintf(fid,'%d %d \n',2,Mobj.nVerts);
  for i=1:Mobj.nVerts
    fprintf(fid,'%f %f \n',Mobj.lon(i),Mobj.lat(i));
  end;
  
  % elements
  fprintf(fid,'Triangles\n');
  fprintf(fid,'%d %d\n',3,Mobj.nElems);
  for i=1:Mobj.nElems
    fprintf(fid,'%d %d %d \n',Mobj.tri(i,1)-1,Mobj.tri(i,2)-1,Mobj.tri(i,3)-1);
  end;
else
  % open datatank file
  fid = fopen(fname,'w');
  fprintf(fid,'DataTank Binary File\n');
  fprintf(fid,'Bathymetry_Cart\n');
  fprintf(fid,'TriangularMesh2D\n');

  %Grid
  fprintf(fid,'Points\n');
  fprintf(fid,'%d %d \n',2,Mobj.nVerts);
  for i=1:Mobj.nVerts
    fprintf(fid,'%f %f \n',Mobj.x(i),Mobj.y(i));
  end;

  % elements
  fprintf(fid,'Triangles\n');
  fprintf(fid,'%d %d\n',3,Mobj.nElems);
  for i=1:Mobj.nElems
    fprintf(fid,'%d %d %d \n',Mobj.tri(i,1)-1,Mobj.tri(i,2)-1,Mobj.tri(i,3)-1);
  end;
end;



% bathymetry 
fprintf(fid,'Values\n');
fprintf(fid,'%d\n',Mobj.nVerts);
for i=1:Mobj.nVerts
  fprintf(fid,'%d \n',Mobj.h(i));
end;


% time step 
fprintf(fid,'Time Step\n');
fprintf(fid,'TriangularMesh2D\n');
fprintf(fid,'Use Last Grid\n');
fprintf(fid,'Values\n');
fprintf(fid,'%d\n',Mobj.nVerts);
for i=1:Mobj.nVerts
  fprintf(fid,'%d \n',Mobj.ts(i));
end;

% open boundary nodes
if(Mobj.nObs > 0)
fprintf(fid,'ObcNodes\n');  
fprintf(fid,'NumberList\n');  
fprintf(fid,'%d\n',sum(Mobj.nObcNodes)); 
for i=1:Mobj.nObcNodes(1) 
  fprintf(fid,'%d\n',Mobj.obc_nodes(1,i)); 
end;
end;

% river nodes
if(Mobj.nRivers > 0)
fprintf(fid,'RivNodes\n');  
fprintf(fid,'NumberList\n');  
fprintf(fid,'%d\n',Mobj.nRivNodes(1)); 
for i=1:Mobj.nRivNodes(1); 
  fprintf(fid,'%d\n',Mobj.riv_nodes(1,i)); 
end;
end;
  
fclose(fid);

