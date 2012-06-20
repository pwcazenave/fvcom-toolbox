function [Mobj]  = setup_metrics(Mobj)

% Setup metrics for mesh object Mesh  
%
% [Mobj] = setup_metrics(Mobj)
%
% DESCRIPTION:
%    Setup metrics for secondary connectivity (nodes surrounding nodes) for Mesh
%
% INPUT
%    Mobj = Matlab mesh object
%
% OUTPUT:
%    Mobj = Matlab mesh object with secondary connectivity
%
% EXAMPLE USAGE
%    Mobj = setup_metrics(Mobj)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
subname = 'setup_metrics';
global ftbverbose
if(ftbverbose)
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;


%------------------------------------------------------------------------------
% Calculate metrics
%------------------------------------------------------------------------------

% set local arrays
tri    = Mobj.tri;
nElems = Mobj.nElems;
nVerts = Mobj.nVerts;

% set xc/yc if they don't exist
if(~isfield(Mobj,'xc'));
  Mobj.xc = zeros(Mobj.nElems,1);
  Mobj.yc = zeros(Mobj.nElems,1);
  for i=1:Mobj.nElems
    Mobj.xc(i) = sum(Mobj.x(Mobj.tri(i,1:3)))/3.;
    Mobj.yc(i) = sum(Mobj.y(Mobj.tri(i,1:3)))/3.;
  end;
end;


% determine edges
nEdges = nElems*3;
edge = zeros(nEdges,2);
icnt = 1;
for i=1:nElems
  edge(icnt  ,1:2) = tri(i,1:2);
  edge(icnt+1,1:2) = tri(i,2:3);
  edge(icnt+2,1:2) = tri(i,[3,1]);
  icnt = icnt + 3;
end;

% determine nodes surrounding nodes (no specific order) 
ntsn = zeros(nVerts,1);
nbsn = zeros(nVerts,12);

for i=1:nEdges
  i1 = edge(i,1);
  i2 = edge(i,2);
  [lmin,loc] = min(abs(nbsn(i1,:)-i2));
  if(lmin ~= 0);
    ntsn(i1) = ntsn(i1)+1;
    nbsn(i1,ntsn(i1)) = i2;
  end;
  [lmin,loc] = min(abs(nbsn(i2,:)-i1));
  if(lmin ~= 0);
    ntsn(i2) = ntsn(i2)+1;
    nbsn(i2,ntsn(i2)) = i1;
  end;
end;



% transfer to struct
Mobj.ntsn = ntsn;
Mobj.nbsn = nbsn;
Mobj.have_mets = true;

if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end;

