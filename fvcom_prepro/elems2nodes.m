function [fieldout]  = elems2nodes(fieldin,Mobj)

% Transfer a field from elements to vertices  
%
% [fieldout] = function elems2nodes(fieldin,Mobj)  
%
% DESCRIPTION:
%    Smooth a vertex based field 
%
% INPUT
%    Mobj         = Matlab mesh object
%    fieldin      = element-based field
%
% OUTPUT:
%    fieldout = vertex-based field
%
% EXAMPLE USAGE
%    f = smoothfield(fc,Mobj)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
subname = 'elems2nodes';
%fprintf('\n')
%fprintf(['begin : ' subname '\n'])

%------------------------------------------------------------------------------
% Parse input
%------------------------------------------------------------------------------

if(exist('fieldin')*exist('Mobj') == 0)
	error('arguments to elems2nodes are missing')
end;

if(length(fieldin) ~= Mobj.nElems)
	error('field size in elems2nodes is not the same as number of elements in Mesh')
end;

%------------------------------------------------------------------------------
% Tranfser
%------------------------------------------------------------------------------
fieldout = zeros(Mobj.nVerts,1);
count    = zeros(Mobj.nVerts,1);

for i=1:Mobj.nElems
	n1 = Mobj.tri(i,1);
	n2 = Mobj.tri(i,2);
	n3 = Mobj.tri(i,3);
	fieldout(n1) = fieldout(n1) + fieldin(i);  count(n1) = count(n1) + 1;
	fieldout(n2) = fieldout(n2) + fieldin(i);  count(n2) = count(n2) + 1;
	fieldout(n3) = fieldout(n3) + fieldin(i);  count(n3) = count(n3) + 1;
end;
fieldout = fieldout./real(count);

%fprintf(['end   : ' subname '\n'])

