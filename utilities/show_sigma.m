%function [] = show_sigma(meshfile,bathfile,sigmafile) 
% plot a sigma distribution along a user-selected line
%
% function [] = show_sigma(meshfile,bathfile,sigmafile) 
%
% DESCRIPTION:
%    plot a sigma distribution along a user-selected line
%
% INPUT:
%    meshfile:   fvcom grid file
%    bathfile:   fvcom bathymetry file
%    sigmafile:  fvcom sigma distribution file
%
% OUTPUT:
%
% EXAMPLE USAGE
%    show_sigma('tst_grd.dat','tst_dep.dat','sigma.dat') 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
meshfile = 'tst_grd.dat';
bathfile = 'tst_dep.dat';
sigmafile = 'sigma.dat';
close all;

% read the mesh 
fid = fopen(meshfile,'r');
if(fid  < 0)
  error(['file: ' meshfile ' does not exist']);
end;
C = textscan(fid, '%s %s %s %d', 1);
Nverts = C{4};
C = textscan(fid, '%s %s %s %d', 1);
Nelems = C{4};
x = zeros(Nverts,3);
tri = zeros(Nelems,3);
fprintf('Nverts: %d\n',Nverts)
fprintf('Nelems: %d\n',Nelems)
for i=1:Nelems 
  C = textscan(fid, '%d %d %d %d %d', 1);
  tri(i,1) = C{2}; 
  tri(i,2) = C{3}; 
  tri(i,3) = C{4}; 
end;
for i=1:Nverts 
  C = textscan(fid, '%d %f %f %f ', 1);
  x(i,1) = C{2}; 
  x(i,2) = C{3}; 
end;

% read the bathy 
fid = fopen(bathfile,'r');
if(fid  < 0)
  error(['file: ' bathfile ' does not exist']);
end;
C = textscan(fid, '%s %s %s %d', 1);
for i=1:Nverts
  C = textscan(fid, '%f %f %f ', 1);
  x(i,3) = C{3};
end;

fprintf('min topography %f\n',min(x(:,3)));
fprintf('max topography %f\n',max(x(:,3)));

% read the sigma file
fid = fopen(sigmafile,'r');
if(fid  < 0)
  error(['file: ' sigmafile ' does not exist']);
end;
C = textscan(fid, '%s %s %s %s %s %d', 1);
nlev = C{6};
C = textscan(fid, '%s %s %s %s %s', 1);
sigtype = char(C{5});
C = textscan(fid, '%s %s %f', 1);
du = C{3};
C = textscan(fid, '%s %s %f', 1);
dl = C{3};
C = textscan(fid, '%s %s %s %s %f', 1);
min_constant_depth = C{5};
C = textscan(fid, '%s %s %d', 1);
ku = C{3};
C = textscan(fid, '%s %s %d', 1);
kl = C{3};
C = textscan(fid, '%s %s %f %f %f %f %f %f %f %f %f %f %f', 1);
for i=1:ku
 zku(i) = C{2+i};
end;
C = textscan(fid, '%s %s %f %f %f %f %f %f %f %f %f %f %f', 1);
for i=1:kl
 zkl(i) = C{2+i};
end;

fprintf('nlev %d\n',nlev)
fprintf('sigtype %s\n',sigtype)
fprintf('du %d\n',du)
fprintf('dl %d\n',dl)
fprintf('min_constant_depth %f\n',min_constant_depth)
fprintf('ku %d\n',ku)
fprintf('kl %d\n',kl)
fprintf('zku %d\n',zku)
fprintf('zkl %d\n',zkl)



% generate the sigma coordinates
fprintf('select two endpoints of a transect')
patch('Vertices',[x(:,1),x(:,2)],'Faces',tri,...
       'Cdata',x(:,3),'edgecolor','interp','facecolor','interp');
axis equal

% plot to get a line
%fprintf('select two end points of a transect with your mouse');
[xt,yt] = ginput(2);
hold on

npts = 25;
ds = (xt(2)-xt(1))/(npts-1);
xline = xt(1):ds:xt(2);
ds = (yt(2)-yt(1))/(npts-1);
yline = yt(1):ds:yt(2);
plot(xline,yline,'w+')
sline(1) = 0.;
for i=2:npts
  sline(i) = sline(i-1) + sqrt(   (xline(i)-xline(i-1))^2 + (yline(i)-yline(i-1))^2);
end;

% interpolate the bathymetry along the line
[zline] = griddata(x(:,1),x(:,2),x(:,3),xline,yline);

figure
plot(sline,-zline)

% generate the sigma coordinates along the line
xslice = zeros(npts,nlev);
yslice = zeros(npts,nlev);

% calculate the sigma distributions along the transect
if(max(sigtype(1:3)=='gen') || max(sigtype(1:3)=='GEN'))
  for i=1:npts
     z(i,1:nlev) = sigma_gen(nlev,dl,du,kl,ku,zkl,zku,zline(i),min_constant_depth);
  end;
elseif(max(sigtype(1:3)=='uni') || max(sigtype(1:3)=='UNI'))
  for i=1:npts
    z(i,1:nlev) = 0:-1/double(nlev-1):-1; 
  end;
else
  error('cant do that sigtype')
end;



for i=1:npts
  xslice(i,1:nlev) = sline(i);
  yslice(i,1:nlev) = z(i,1:nlev)*zline(i);
end;


% plot the mesh along the transect
for k=1:nlev
  plot(xslice(:,k),yslice(:,k),'k'); hold on;
end;
for k=1:npts
  plot(xslice(k,:),yslice(k,:),'k'); hold on;
end;
axis([xslice(1,1),xslice(end,1),min(yslice(:,end)),5])
title('sigma distribution along the transect');


