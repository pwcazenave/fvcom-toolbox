function [tri,x,y] = tri_from_bndry(xb,yb,path_to_triangle); 
% tesselate a triangular mesh using "triangle" from a single, unclosed boundary 
%
% input:  xb,yb (boundary nodes, unclosed)
%         path_to_triangle:  full path to triangle executable ( a string )
% output: tri,x,y (connectivity, nodes)


% dump poly file
fid = fopen('junk.poly','w');
npts = length(xb);
fprintf(fid,'%d %d %d %d\n',npts,2,1,1);
for i=1:npts
  fprintf(fid,'%d %f %f %d %d\n',i,xb(i),yb(i),1,1);
end;
fprintf(fid,'%d %d\n',npts,1);
for i=1:npts-1
  fprintf(fid,'%d %d %d %d\n',i,i,i+1,1);  
end;
fprintf(fid,'%d %d %d %d\n',npts,npts,1,1);  
fprintf(fid,'%d\n',0);
fclose(fid);

% triangulate
system([path_to_triangle ' -pqm -Q junk.poly']);

% read in the triangulate file 
fid = fopen('junk.1.node');
vec = fgetl(fid);
veci = sscanf(vec,'%d %d %d %d');
nvts = veci(1);
b = fscanf(fid,'%d %f %f %d %d',nvts*5); 
a = reshape(b,5,nvts);
x = a(2,:);
y = a(3,:);
fclose(fid);

fid = fopen('junk.1.ele');
vec = fgetl(fid);
veci = sscanf(vec,'%d %d %d');
nelems = veci(1);
b = fscanf(fid,'%d %d %d %d',nelems*4); 
a = reshape(b,4,nelems);
tri = zeros(nelems,3); 
tri(:,1) = a(2,:); 
tri(:,2) = a(3,:); 
tri(:,3) = a(4,:); 
