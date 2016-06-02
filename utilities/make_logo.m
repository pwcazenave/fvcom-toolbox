% Make the toolbox logo. Adapted from:
%   mathworks.com/help/matlab/examples/creating-the-matlab-logo.html
%
% Note: this will not produce the same figure each time as there's a call
% to rand in the script.
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory).
% 
% ChangeLog:
%   2016-06-02 First version.

ftbverbose = 1;

[base, subname] = fileparts(mfilename('fullpath'));
cd(base)

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

% Make an initial membrane and then triangulate it and plot that.

n = 5; % size of the membrane
L = 160 * membrane(1, n);
[X, Y] = meshgrid(1:n*2 + 1, 1:n*2 + 1);

% Interpolate onto an irregular grid.
XX = X .* (0.95 + mod(abs(rand(size(X))), 0.05));
YY = Y .* (0.95 + mod(abs(rand(size(Y))), 0.05));
Z = griddata(X, Y, L, XX, YY);

tri = delaunayTriangulation(XX(:), YY(:));

% Do the plot replicating (as closely as possible) the MATLAB logo.
close all

f = figure;
ax = axes;
s = trisurf(tri.ConnectivityList, XX(:), YY(:), Z(:));
s.EdgeColor = 'none';
view(3)
ax.XLim = [1 max(X(:))];
ax.YLim = [1 max(Y(:))];
ax.ZLim = [min(Z(:)), max(Z(:))];
ax.CameraUpVector = [0 0 1];
ax.CameraViewAngle = 8;
ax.Position = [0 0 1 1];
ax.DataAspectRatio = [1 1 20];
l1 = light;
l1.Position = [160 400 80];
l1.Style = 'local';
l2 = light;
l2.Position = [.5 -1 .4];
l2.Color = [0.8 0.8 0];
caxis([min(Z(:)), max(Z(:))])
s.FaceLighting = 'gouraud';
s.AmbientStrength = 0.3;
s.DiffuseStrength = 0.9;
s.BackFaceLighting = 'lit';
s.SpecularStrength = 1;
s.SpecularColorReflectance = 1;
s.SpecularExponent = 7;
axis off
colormap(hot)

% Transparent background and save to PDF.
set(gca,'color','none') 
print(gcf, '-dpdf', '../doc/fvcom-toolbox.pdf', '-painters')

if ftbverbose
    fprintf('end   : %s \n', subname)
end
