function [spongeRadius] = calc_sponge_radius(Mobj,Nlist) 

% Calculate a variable sponge radius based on distance to the boundary
% node's furthest neighbour.
% (Adapted from Phil Hall's 'produce_netcdf_input_data.py')
%
% spongeRadius = calc_sponge_radius(Mobj,Nlist) 
%
% DESCRIPTION
%    Calculates the sponge radius for each node on the open boundary, based
%    on the minimum of either the distance to the node's furthest
%    neighbour, or 100 km.
%
% INPUT
%    Mobj = Matlab mesh object
%    Nlist = List of nodes
%
% OUTPUT
%    spongeRadius = List of variable sponge radii
%
% EXAMPLE USAGE
%    spongeRadius = calc_sponge_radius(Mobj,Nlist)
%
% Author(s)
%    Karen Amoudry (National Oceanography Centre, Liverpool)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%    2013-01-02 KJA bug fix: amended usage of 'unique' to prevent it from
%    sorting the values it returns. Amended by Pierre to support pre-2012
%    versions of MATLAB whilst giving the same result.
%    2013-07-16 KJA: adapted function to remove dependency on the Matlab
%    Mapping Toolbox.
%
%==========================================================================
subname = 'calc_sponge_radius';
global ftbverbose
if(ftbverbose)
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end

%--------------------------------------------------------------------------
% Get a unique list and make sure they are in the range of node numbers 
%--------------------------------------------------------------------------
% Make this work in versions of MATLAB older than 2012a (newer versions
% can just use unique(A, 'stable'), but checking versions is a pain).
[~, Nidx] = unique(Nlist);
Nlist = Nlist(sort(Nidx));

spongeRadius = 100000+zeros(size(Nlist));

% For each node on the open boundary
for i =1:length(Nlist)
    % Find the neighbouring nodes
    [r,c]=find(Mobj.tri==Nlist(i));
    neighbours = Mobj.tri(r,:);
    [~,neighidx] = unique(Mobj.tri(r,:));
    neighbours = neighbours(sort(neighidx));
    
    % Remove the node of interest from the neighbours list
    n = find(neighbours~=Nlist(i));
    neighbours = neighbours(n);
    
    % Calculate the arc length (in degrees) between the node and its
    % neighbours
%     arclen = distance(Mobj.lat(Nlist(i)),Mobj.lon(Nlist(i)),...
%         Mobj.lat(neighbours),Mobj.lon(neighbours));
%     % Convert from degrees to whole metres
%     arclen = ceil(1000*deg2km(arclen));
    
    
    % Adapted to avoid using Mapping toolbox
    % Step 1: convert lat/lon to radians
    lat1 = Mobj.lat(Nlist(i)) .* pi./180;
    lon1 = Mobj.lon(Nlist(i)) .* pi./180;
    lat2 = Mobj.lat(neighbours) .* pi./180;
    lon2 = Mobj.lon(neighbours) .* pi./180;
    % Step 2: calculate distance in radians
    a = sin((lat2-lat1)/2).^2 + cos(lat1) .* cos(lat2) .*...
        sin((lon2-lon1)/2).^2;
    arclen = 2 * atan2(sqrt(a),sqrt(1 - a));

    % Calculate distance in whole metres
    arclen = ceil(1000 .* 6371 .* arclen);
    
    % If the smallest distance is less than 100km, keep it
    if min(arclen)<spongeRadius(i)
        spongeRadius(i)=min(arclen);
    end
end

if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end

