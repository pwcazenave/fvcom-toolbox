function [Mobj]  = add_stations_list(Mobj,Positions,Names,Dist)

% Add a set of stations at which FVCOM will output time series.
%
% [Mobj] = add_stations_list(Mobj,Positions,Names,Dist)
%
% DESCRIPTION:
%    Supply a list of positions (in the same coordinate system as the
%    native coordinates of the grid) and a cell array of names. Nearest
%    grid node to those supplied will be used in the output file.
%
% INPUT
%    Mobj = Matlab mesh object
%    Positions = 2xn array of the XY positions of the stations
%    Names = Cell array of the names of the stations defined in Positions
%    Dist = Maximum distance from a station for a node to be included
%
%    Optionally supply positions as a 4xn array with spherical x and y and
%    cartesian x and y in columns 1, 2, 3 and 4, respectively. The
%    values in Mobj.nativecoords will be used for the distance check, so
%    ensure Dist is in those units.
%
% OUTPUT:
%    Mobj = Matlab mesh object with an additional cell array (stations)
%    containing:
%       {id, x, y, nodelist, depth, 'station name', elem}
%    where id is the station ID (from 1 number of valid stations, x and y
%    are the coordinates of the station, nodelist is the indices of the
%    model grid nodes, depth is the depth at the grid node, station name is
%    a string of the name from the input file and eleme is the grid element
%    ID closest to each station.
%
% EXAMPLE USAGE
%    Mobj = add_stations_list(Mobj, [-5.54, 50.103; -3.0865, 58.441], ...
%    {'Newlyn', 'Wick'}, 0.25)
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2012-11-30 First version.
%    2015-01-14 Add support for exporting the element ID closest to the
%    station of interest (as well as the node ID which was already there).
%    Fix some formatting issues.
%
%==========================================================================
subname = 'add_stations_list';
global ftbverbose
if ftbverbose
  fprintf('\nbegin : %s\n', subname)
end

%--------------------------------------------------------------------------
% Check the inputs
%--------------------------------------------------------------------------
nPos = size(Positions, 1);
nNames = size(Names, 1);
if nPos ~= nNames
    error('The number of the supplied station positions and names do not match (%i and %i respectively)', nPos, nNames)
end

%--------------------------------------------------------------------------
% For each site in the supplied positions, find the nearest node ID
%--------------------------------------------------------------------------

% Check for whether the input has both spherical and cartesian.
if size(Positions, 2) > 2
    % Now check for which is the native coordinate system, and output the
    % station positions in that coordinate system.
    if strcmpi(Mobj.nativeCoords, 'cartesian')
        cols = [3, 4];
    elseif strcmpi(Mobj.nativeCoords, 'spherical')
        cols = [1, 2];
    else
        error('Unknown native coordinate system string: %s', Mobj.nativeCoords)
    end
else
    % We have to assume the positions are in the grid's native coordinate
    % system.
    cols = [1, 2];
end

inc = 1;
out = cell(1); % don't preallocate as we don't know how many we'll have

for s = 1:nPos
    [node, dist] = find_nearest_pt(Positions(s, cols(1)), Positions(s, cols(2)), Mobj);
    [~, elem] = min(abs(sqrt((Mobj.xc - Positions(s, cols(1))).^2 + Mobj.yc - Positions(s, cols(2))).^2));

    if dist >= Dist
        % Skip out for this station
        if ftbverbose
            fprintf('Skipping station %s (%g, %g). Nodal distance from station position = %f\n', Names{s}, Positions(s, 1), Positions(s, 2), dist)
        end
        continue
    end
    out{inc} = {inc, Positions(s, cols(1)), Positions(s, cols(2)), node, Mobj.h(node), Names{s}, elem};
    inc = inc + 1;
end

if ~isempty(out)
    Mobj.stations = out;
else
    Mobj.stations = [];
    if ftbverbose
        fprintf('No stations found within the model domain.\n')
    end
end

if ftbverbose
  fprintf('\nend   : %s\n', subname)
end