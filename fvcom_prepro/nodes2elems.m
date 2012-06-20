function [fieldout]  = nodes2elems(fieldin,Mobj)

% Transfer a field from vertices to elements
%
% function [fieldout] = nodes2elems(fieldin,Mobj)  
%
% DESCRIPTION:
%    Smooth a vertex based field 
%
% INPUT
%    Mobj         = Matlab mesh object
%    fieldin      = vertex-based field
%
% OUTPUT:
%    fieldout = element-based field
%
% EXAMPLE USAGE
%    f = smoothfield(fv,Mobj)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
subname = 'nodes2elems';
global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

%------------------------------------------------------------------------------
% Parse input
%------------------------------------------------------------------------------

if(exist('fieldin')*exist('Mobj') == 0)
	error('arguments to nodes2elems are missing')
end;

if(length(fieldin) ~= Mobj.nVerts)
	error('field size in nodes2elems is not the same as number of nodes in Mesh')
end;

%------------------------------------------------------------------------------
% Tranfser
%------------------------------------------------------------------------------
fieldout = zeros(Mobj.nElems,1);

for i=1:Mobj.nElems
	fieldout(i) = sum(fieldin(Mobj.tri(i,1:3)))/3.; 
end;


if(ftbverbose);
  fprintf(['end   : ' subname '\n'])
end;

