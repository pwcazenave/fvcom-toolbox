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
%    Mobj = Matlab mesh object with an additional cell array containing id,
%    x, y, nodelist, depth and station name.
%
% EXAMPLE USAGE
%    Mobj = add_stations_list(Mobj, [-5.54, 50.103; -3.0865, 58.441], ...
%    {'Newlyn', 'Wick'}, 0.25)
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
%
% Revision history
%    2012-11-30 First version.
%
%==========================================================================
subname = 'add_stations_list';
global ftbverbose
if(ftbverbose)
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

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
% out = cell(nPos, 1);

for s=1:nPos
    [node, dist] = find_nearest_pt(Positions(s, cols(1)), Positions(s, cols(2)), Mobj);

    if dist >= Dist
        % Skip out for this station
        if(ftbverbose)
            fprintf('Skipping station %s (%g, %g). Nodal distance from station position = %f\n', Names{s}, Positions(s, 1), Positions(s, 2), dist)
        end
        continue
    end
    out{inc} = {inc, Positions(s, cols(1)), Positions(s, cols(2)), node, Mobj.h(node), Names{s}};
    inc = inc + 1;
end

Mobj.stations = out;
