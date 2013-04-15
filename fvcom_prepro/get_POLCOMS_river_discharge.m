function Mobj = get_POLCOMS_river_discharge(Mobj, polcoms_flow)
% Parse the POLCOMS rivers data file.
%
% get_POLCOMS_rivers(Mobj, polcoms_file, polcoms_grid, polcoms_ij)
%
% DESCRIPTION:
%   Takes POLCOMS grid and index files and returns the positions and names
%   of the rivers within the index file.
%
% INPUT:
%   Mobj - MATLAB mesh object into which the outputs will be added.
%   polcoms_flow - flow data file(s) from POLCOMS. For multiple files, give
%       in chronological order as a cell array of file names.
% 
% OUTPUT:
%   Mobj.river_discharge - array of discharges for the rivers in the
%   POLCOMS river flow file(s).
%
% EXAMPLE USAGE:
%   Mobj = get_POLCOMS_river_discharge(Mobj, 'polcoms.flw')
% 
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2013-04-15 - First version from the original get_POLCOMS_rivers.m.
%
%==========================================================================

subname = 'get_POLCOMS_river_discharge';

global ftbverbose;
if ftbverbose
    fprintf(['\nbegin : ' subname '\n'])
end

% Check inputs
if exist(polcoms_flow, 'file') ~= 2
    error('file: %s does not exist or is not readable.', polcoms_grid)
end

% The POLCOMS river file has a pretty straightforward format of a 2D array
% of river along x and time along y. Since it's a simple text file with no
% weird format, we'll just read it in with load.
if iscell(polcoms_flow)
    pc_riv = [];
    for rr = 1:length(polcoms_flow)
        pc_riv = [pc_riv; load(polcoms_flow{rr})];
    end
    clear rr
else
    pc_riv = load(polcoms_flow);
end
[pc_nt, pc_nr] = size(pc_riv);

% Return the number of river and the number of times along with the actual
% discharge data.
Mobj.rivers.num = pc_nr;
Mobj.rivers.num_time = pc_nt;
Mobj.rivers.discharge = pc_riv;

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end
