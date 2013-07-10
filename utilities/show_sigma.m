function show_sigma(meshfile,bathfile,sigmafile,varargin)
% plot a sigma distribution along a user-selected line
%
% show_sigma(meshfile,bathfile,sigmafile) 
%
% DESCRIPTION:
%    plot a sigma distribution along a user-selected line
%
% INPUT:
%    meshfile:   fvcom grid file
%    bathfile:   fvcom bathymetry file
%    sigmafile:  fvcom sigma distribution file
%    npts:       number of points in the user-selected line (optional,
%                defaults to 25).
%
% OUTPUT:
%
% EXAMPLE USAGE
%    show_sigma('tst_grd.dat', 'tst_dep.dat', 'sigma.dat', 50)
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2013-01-08 Added more resilient reading of the sigma coordinates file
%   (can now handle comments and blank lines). Also add extra (optional)
%   argument to increase the profile sampling frequency.
%   
%==============================================================================

close all

global ftbverbose

if nargin == 4
    npts = varargin{1};
else
    npts = 25;
end

% read the mesh
fid = fopen(meshfile,'r');
if(fid  < 0)
    error(['file: ' meshfile ' does not exist']);
end
C = textscan(fid, '%s %s %s %d', 1);
Nverts = C{4};
C = textscan(fid, '%s %s %s %d', 1);
Nelems = C{4};
x = zeros(Nverts,3);
tri = zeros(Nelems,3);
if ftbverbose
    fprintf('Nverts: %d\n',Nverts)
    fprintf('Nelems: %d\n',Nelems)
end

for i=1:Nelems
    C = textscan(fid, '%d %d %d %d %d', 1);
    tri(i,1) = C{2};
    tri(i,2) = C{3};
    tri(i,3) = C{4};
end
for i=1:Nverts
    C = textscan(fid, '%d %f %f %f ', 1);
    x(i,1) = C{2};
    x(i,2) = C{3};
end

% read the bathy
fid = fopen(bathfile,'r');
if(fid  < 0)
    error(['file: ' bathfile ' does not exist']);
end
C = textscan(fid, '%s %s %s %d', 1);
for i=1:Nverts
    C = textscan(fid, '%f %f %f ', 1);
    x(i,3) = C{3};
end

if ftbverbose
    fprintf('min topography %f\n',min(x(:,3)));
    fprintf('max topography %f\n',max(x(:,3)));
end

% read the sigma file
fid = fopen(sigmafile,'r');
if(fid  < 0)
    error(['file: ' sigmafile ' does not exist']);
end

while ~feof(fid)
    line = fgetl(fid);
    if isempty(line) || strncmp(line, '!', 1) || ~ischar(line)
        continue
    end
    key = lower(line(1:3));
    C = strtrim(regexpi(line, '=', 'split'));
    switch key
        case 'num'
            nlev = str2double(C{2});
        case 'sig'
            sigtype = C{2};
        case 'du '
            du = str2double(C{2});
        case 'dl '
            dl = str2double(C{2});
        case 'min'
            min_constant_depth = str2double(C{2});
        case 'ku '
            ku = str2double(C{2});
        case 'kl '
            kl = str2double(C{2});
        case 'zku'
            s = str2double(regexp(C{2}, ' ', 'split'));
            zku = zeros(ku, 1);
            for i = 1:ku
                zku(i) = s(i);
            end
        case 'zkl'
            s = str2double(regexp(C{2}, ' ', 'split'));
            zkl = zeros(kl, 1);
            for i = 1:kl
                zkl(i) = s(i);
            end
    end
end

if ftbverbose
    fprintf('nlev %d\n',nlev)
    fprintf('sigtype %s\n',sigtype)
    fprintf('du %d\n',du)
    fprintf('dl %d\n',dl)
    fprintf('min_constant_depth %f\n',min_constant_depth)
    fprintf('ku %d\n',ku)
    fprintf('kl %d\n',kl)
    fprintf('zku %d\n',zku)
    fprintf('zkl %d\n',zkl)
end

% generate the sigma coordinates

% fix "java.lang.IllegalArgumentException: adding a container to a container
% on a different GraphicsDevice" error. See
% http://www.mathworks.co.uk/matlabcentral/newsreader/view_thread/169024.
set(0,'DefaultFigureRenderer','OpenGL')

figure(1)
patch('Vertices',[x(:,1),x(:,2)],'Faces',tri,...
    'Cdata',x(:,3),'edgecolor','interp','facecolor','interp');
axis equal

% plot to get a line
fprintf('select two end points of a transect with your mouse... ');
[xt,yt] = ginput(2);
hold on
fprintf('done.\n')

ds = (xt(2)-xt(1))/(npts-1);
xline = xt(1):ds:xt(2);
ds = (yt(2)-yt(1))/(npts-1);
yline = yt(1):ds:yt(2);
plot(xline, yline, 'w+')
sline = zeros(1, npts);
for i=2:npts
    sline(i) = sline(i-1) + sqrt((xline(i)-xline(i-1))^2 + (yline(i)-yline(i-1))^2);
end

% interpolate the bathymetry along the line
zline = griddata(x(:,1), x(:,2), x(:,3), xline, yline);

figure(2)
plot(sline, -zline)

% generate the sigma coordinates along the line
xslice = zeros(npts,nlev);
yslice = zeros(npts,nlev);
z = nan(npts, nlev);

% calculate the sigma distributions along the transect
switch lower(sigtype)
    case 'generalized'
        for i=1:npts
            z(i,1:nlev) = sigma_gen(nlev,dl,du,kl,ku,zkl,zku,zline(i),min_constant_depth);
        end
    case 'uniform'
        for i=1:npts
            z(i,1:nlev) = 0:-1/double(nlev-1):-1;
        end
    otherwise
        error('Can''t do that sigtype')
end

for i=1:npts
    xslice(i,1:nlev) = sline(i);
    yslice(i,1:nlev) = z(i,1:nlev) * zline(i);
end


% plot the mesh along the transect
for k=1:nlev
    plot(xslice(:,k),yslice(:,k),'k'); hold on;
end
for k=1:npts
    plot(xslice(k,:),yslice(k,:),'k'); hold on;
end
axis([xslice(1,1),xslice(end,1),min(yslice(:,end)),5])
title('sigma distribution along the transect');
