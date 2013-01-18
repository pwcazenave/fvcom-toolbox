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
%    Karen Thurston (National Oceanography Centre, Liverpool)
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
Nlist = unique(Nlist);

spongeRadius = 100000+zeros(size(Nlist));

% For each node on the open boundary
for i =1:length(Nlist)
    % Find the neighbouring nodes
    [r,c]=find(Mobj.tri==Nlist(i));
    neighbours = unique(Mobj.tri(r,:));
    
    % Remove the node of interest from the neighbours list
    n = find(neighbours~=Nlist(i));
    neighbours = neighbours(n);
    
    % Calculate the arc length (in degrees) between the node and its
    % neighbours
    arclen = distance(Mobj.lat(Nlist(i)),Mobj.lon(Nlist(i)),...
        Mobj.lat(neighbours),Mobj.lon(neighbours));
    
    % Convert from degrees to whole metres
    arclen = ceil(1000*deg2km(arclen));
    
    % If the smallest distance is less than 100km, keep it
    if min(arclen)<spongeRadius(i)
        spongeRadius(i)=min(arclen);
    end
end

if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end

