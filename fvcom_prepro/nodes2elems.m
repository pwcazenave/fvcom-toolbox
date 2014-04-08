function fieldout = nodes2elems(fieldin,Mobj)
% Transfer a field from vertices to elements
%
% function fieldout = nodes2elems(fieldin, Mobj)  
%
% DESCRIPTION:
%    Transfer a field from vertices (nodes) to elements
%
% INPUT
%    Mobj         = Matlab mesh object
%    fieldin      = vertex-based field
%
% OUTPUT:
%    fieldout = element-based field
%
% EXAMPLE USAGE
%    f = smoothfield(fv, Mobj)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2014-03-10 Made the calculation array-based instead of loop-based
%   (increases performance).
%   
%==========================================================================

subname = 'nodes2elems';

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

%--------------------------------------------------------------------------
% Parse input
%--------------------------------------------------------------------------

if exist('fieldin', 'var') ~= 1 || exist('Mobj', 'var') ~= 1
	error('arguments to nodes2elems are missing')
end

if length(fieldin) ~= Mobj.nVerts
	error('field size in nodes2elems is not the same as number of nodes in Mesh')
end

%--------------------------------------------------------------------------
% Tranfser
%--------------------------------------------------------------------------

fieldout = mean(fieldin(Mobj.tri), 2);

if ftbverbose
    fprintf('end   : %s \n', subname)
end
