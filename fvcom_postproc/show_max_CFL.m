function [ maxCFLs, fighandle ] = show_max_CFL(grdFile, depFile, ncFile, extTS, coordtype, fig_flag)
%SHOW_MAX-CFL Function to find the max CFL encountered in each mesh element during an FVCOM model
%run. NB need enough RAM to load all u, v and zeta data from the given
%ncFile.
%   Inputs: grdFile, depFile: casename_grd.dat and casename_dep.dat files
%                              from an FVCOM model
%           ncFile : Output from FVCOM. Must include u, v, and zeta.
%           extTS : External timestep at which to calculate CFL, in
%                   seconds.
%           coordtype: 'lonlat' or 'cartesian'. If cartesian, coordinates
%                  should be in metres. If lonlat, they should be in
%                  degrees.
%           fig_flag: Should a figure be plotted? true or false.
%           

% The formulation of CFL used here is taken from the MIKE 3 manual

% Simon Waldman (Marine Scotland Science / Heriot-Watt University), March 2018

%% Check inputs

assert( nargin == 6, 'Wrong number of arguments.');
assert( exist(grdFile, 'file')==2, 'Cannot find file %s', grdFile );
assert( exist(depFile, 'file')==2, 'Cannot find file %s', depFile );
assert( exist(ncFile, 'file')==2, 'Cannot find file %s', ncFile );
assert( ischar(coordtype), 'Coord type must be ''lonlat'' or ''cartesian''.' );
assert( islogical( fig_flag ), 'fig_flag should be logical.' );

switch coordtype
    case 'lonlat'
        lonlat=true;
    case 'latlon'
        lonlat=true;
    case 'cartesian'
        lonlat=false;
    otherwise
        error('Coord type must be ''lonlat'' or ''cartesian''.');
end

% does the ncFile have velocities & water levels?

info = ncinfo(ncFile);
assert( any(strcmp('u', {info.Variables.Name})) && any(strcmp('v', {info.Variables.Name})) && ...
    any(strcmp('zeta', {info.Variables.Name})), 'netCDF file must include the fields u, v and zeta.' );

%% Do the stuff.

disp('Reading mesh & bathymetry...');
M = read_fvcom_mesh( grdFile ); %NB this function puts lon and lat in the x and y fields.
M.h = read_fvcom_bath( depFile );
% calculate element depths from mean of nodes
M.hc = mean(M.h(M.tri),2);

disp('Reading U velocities...');
U = ncread(ncFile, 'u');
disp('Reading V velocities...');
V = ncread(ncFile, 'v');
disp('Reading surface elevations...');
Z = ncread(ncFile, 'zeta');

g = 9.81;

NumTS = size(U, 3);
NumEl = size(U, 1);

% Find water depths for each element at each TS. First calculate this for
% nodes, then average them for elements.
disp('Calculating depths...');
NodeDepths = repmat(M.h, 1, NumTS) + Z; %repmat because Z has a time dimension, M.h doesn't.
for n = 1:3
    tmp(:,:,n) = NodeDepths(M.tri(:,n), :);
end
ElDepths = mean(tmp, 3); %this has dimensions of element x TS.
clear Z NodeDepths tmp; %save some memory

% Find minimum characteristic length for each element. This is approximated by the
% minimum edge length. It could be shorter with really long thin triangles, 
% but the 30 degree internal angle minimum for FVCOM means we shouldn't have those.

disp('Calculating triangle sizes...');
CharLen = nan(NumEl, 1);
for e = 1:NumEl %for each element
    xv = M.x(M.tri(e,:)); %vertices.
    yv = M.y(M.tri(e,:));
    %close the triangle by copying the first vertex to the end
    xv(4) = xv(1);
    yv(4) = yv(1);
    %find the edge lengths
    if lonlat
        for a = 1:3
            edges(a) = haversine(yv(a), xv(a), yv(a+1), xv(a+1));
        end
    else
        edges = sqrt( diff(xv).^2 + diff(yv).^2 );
    end
    CharLen(e) = min(edges);
end

% For each element and TS, find the layer with the highest U and V magnitudes.
% Technically doing this with each component rather than the vector magnitude 
% is wrong, but it'll usually be close and where wrong, it'll overestimate the CFL, so it's safe.
maxU = squeeze(max(abs(U), [], 2));
maxV = squeeze(max(abs(V), [], 2));
    
% Find CFL for each element at each TS
CFL = ( 2 .* sqrt( g .* ElDepths ) + maxU + maxV ) .* repmat( (extTS ./ CharLen), 1, NumTS );
%This is based on equation 6.1 on pg 33 of the MIKE hydrodynamic module
%manual (modified for using a single characteristic length rather than
%deltaX/deltaY)

% find the max over time for each element
MaxCFL = max(CFL, [], 2);

[val, I] = max(MaxCFL);
fprintf('Max CFL reached with an external timestep of %.2f secs is approx. %.3f, in Element %i.\n', extTS, val, I);

% find how long the timestep probably could go. We set the CFL to 0.8 and
% apply the same formula backwards.
TargetCFL = 0.8;
MaxTSs = repmat( (TargetCFL .* CharLen), 1, NumTS ) ./ ( 2 .* sqrt( g .* ElDepths ) + maxU + maxV );
%that's still per element per TS. We care about what the smallest is.
OverallMaxTS = min(min(MaxTSs));

fprintf('Max external timestep to reach CFL of %.1f with this mesh would be approx. %.2f seconds.\n', TargetCFL, OverallMaxTS );

% Optionally, plot figure of this
if fig_flag
    CFLfig = figure;
    p = patch();
    p.Vertices = [M.x M.y];
    p.Faces = M.tri;
    p.CData = MaxCFL;
    p.FaceColor = 'flat';
    %p.EdgeColor = [0 0 0];
    p.EdgeColor = 'none';
    p.EdgeAlpha = 0.1;
    p.LineWidth = 0.1;
    cb = colorbar;
    colormap(parula);
    caxis([0 max(MaxCFL)]);
    ylabel(cb, 'Max CFL encountered.');
    axis equal;
    
    % add to the plot markers for the worst n elements
    n = 10; %number of cells to highlight
    [ ~, WorstEls ] = sort(MaxCFL, 'descend');
    WorstElsX = mean(M.x(M.tri(WorstEls(1:n),:)), 2); %finding element centroid coords.
    WorstElsY = mean(M.y(M.tri(WorstEls(1:n),:)), 2);
    hold on;
    plot(WorstElsX, WorstElsY, 'or');
    fprintf('Red circles on plot show the %i mesh elements with the highest CFL.\n', n);
end

%return values
maxCFLs = MaxCFL;
fighandle = CFLfig;

end



function [distm]=haversine(lat1,lon1,lat2,lon2)
% Haversine function to calculate first order distance measurement. Assumes
% spherical Earth surface. Lifted from:
%
% http://www.mathworks.com/matlabcentral/fileexchange/27785

% Could be done more accurately with Mapping Toolbox tools, but don't want
% to require that, and we don't need amazing accuracy.

lat1 = deg2rad(lat1);
lat2 = deg2rad(lat2);
lon1 = deg2rad(lon1);
lon2 = deg2rad(lon2);
R = 6371000;                    % Earth's mean radius in metres
delta_lat = lat2 - lat1;        % difference in latitude
delta_lon = lon2 - lon1;        % difference in longitude
a = sin(delta_lat/2)^2 + cos(lat1) * cos(lat2) * ...
    sin(delta_lon/2)^2;
c = 2 * atan2(sqrt(a), sqrt(1-a));
distm = R * c;                     % distance in metres
end