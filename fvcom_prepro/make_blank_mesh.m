function [Mobj] = make_blank_mesh

% Make a blank mesh object with default params  
%
% [Mobj] = function make_blank_mesh
%
% DESCRIPTION:
%    Make a blank Matlab mesh object
%
% INPUT:
%
% OUTPUT:
%    Mobj = matlab structure containing default
%
% EXAMPLE USAGE
%    Mobj = make_blank_mesh()
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
global ftbverbose
subname = 'make_blank_mesh';
if(ftbverbose)
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

%------------------------------------------------------------------------------
% Set defaults
%------------------------------------------------------------------------------


% dimensions
Mobj.nVerts  = 0;
Mobj.nElems  = 0;
Mobj.nRivers = 0;
Mobj.nObs    = 0;
Mobj.nSponge = 0;
Mobj.riv_nodes = zeros(50,10);
Mobj.obc_nodes = zeros(10,500);
Mobj.sponge_nodes = zeros(10,500);

% flags
Mobj.nativeCoords = 'cartesian';
Mobj.have_lonlat  = false;
Mobj.have_cor     = false;
Mobj.have_xy      = false;
Mobj.have_bath    = false;
Mobj.have_mets    = false;


if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end;


