function RX1matrix=ComputeMatrixRx1_nodes(Z_w, Mobj)
% Computes rx1 matrix as in ComputeMatrixRx1_V2.m downloaded from
% https://github.com/dcherian/tools
%
% function [RX1matrix] = ComputeMatrixRx1_nodes(Z_w, Mobj)
%
% DESCRIPTION:
%    Calculates the hydrostatic consistency condition:
%      r = abs(Sigma/H deltaxH/deltaSigma)
%    this reflects the errors associated with horizontal pressure gradients
%    calculations that are associated with steep bathyemtry and low
%    vertical resolution.
%
% INPUT
%   Z_w    = This is the sigmal layer vertical distribution in Z coordinates
%   Mobj   = needs triangulation and mesh information
%            table. read_sms_mesh provides everything it needs.
%
%
% OUTPUT:
%    RX1matrix = node based field with values of max(rx1)
%
% EXAMPLE USAGE
%    [RX1matrix] = ComputeMatrixRx1_nodes(Z_w, Mobj)
%
% Author(s):
%    Ricardo Torres (Plymouth Marine Laboratory) based on ComputeMatrixRx1_V2
%
% Revision history
%
%   2018-03-22 First version.
%
%==============================================================================

subname = 'ComputeMatrixRx1_nodes';
global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end
%--------------------------------------------------------------------------
% Estimate dimensins of mesh based on array's sizes
%--------------------------------------------------------------------------

N=size(Z_w,2)-1;
nNodes = size(Z_w,1);
RX1matrix=zeros(nNodes,1);
TR = triangulation(Mobj.tri, [Mobj.x, Mobj.y]);
% loop through all Nodes
for iXi=1:nNodes
    rx1=0;
% Find neighbouring nodes for the node under consideration
    ti = cell2mat(vertexAttachments(TR,iXi));
    vertices = setdiff(unique(TR.ConnectivityList(ti,:)),iXi);
% calculate max rx1 values in the vertical and horizontal among all neighbouring nodes
    for nn=1:length(vertices)
        for i=1:N
            a1=abs(Z_w(iXi,i+1) - Z_w(vertices(nn),i+1) + ...
                Z_w(iXi,i) - Z_w(vertices(nn),i));
            b1=abs(Z_w(iXi,i+1) + Z_w(vertices(nn),i+1) - ...
                Z_w(iXi,i) - Z_w(vertices(nn),i));
            quot=abs(a1/b1);
            rx1=max(rx1, quot);
        end
    end
    RX1matrix(iXi)=rx1;
end
if ftbverbose
  fprintf('end   : %s\n', subname)
end
