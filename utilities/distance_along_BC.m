function Mobj=distance_along_BC(Mobj,BCnodes,conf)
% Calculates the distance from coast along the open boundary nodes
%
% Mobj=distance_along_BC(Mobj,BCnodes,conf)
%
% DESCRIPTION:
%    Calculates distance from coast along open boundary nodes and
%    store them in a matlab mesh object
%
% INPUT 
%   Mobj                   = Mesh object structure variable
%   BCnodes                = indices of nodes located at the boundary
%   conf                   = configuration structure variable with the
%   directory where the HJB_Solver_Package is installed
%
% OUTPUT:
%    Mobj = matlab structure containing distance data
%
% EXAMPLE USAGE
%     Mobj=distance_along_BC(Mobj,BCnodes,conf)
% This function needs the HJB_solver package by Shawn Walker and can be downloaded from Matlab central 
% http://www.mathworks.com/matlabcentral/fileexchange/24827-hamilton-jacobi-solver-on-unstructured-triangular-grids

% Author(s):
%    Ricardo Torres  (Plymouth Marine Laboratory)
%
% Revision history
%
%   2015-11-20 First version 
%
%==============================================================================

dump = dbstack;
subname = dump.name;
clear dump
global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end
CD=pwd;
% setup HPJ solver to calculate the distance function for the SMS mesh
 [~,~,~,bnd] = connectivity([Mobj.x,Mobj.y],Mobj.tri);
% remove nodestring from coast.
% BCnodes=[Mobj.read_obc_nodes{:}];
coast_ind=(BCnodes);


 % % calculate distance function 
myParam.Max_Tri_per_Star = 20;
myParam.NumGaussSeidel = 40;
myParam.INF_VAL = 1000000;
myParam.TOL = 1e-12;
% in this case, we will assume the standard Euclidean metric
myMetric = [];
myTM.Vtx=[Mobj.x(:),Mobj.y(:)];
myTM.DoFmap=[Mobj.tri];
myTM.NegMask=false(size(bnd));
myBdy.Nodes=coast_ind(:);
myBdy.Data=zeros(size(myBdy.Nodes));
% 
cd (conf.HJB_Solver_Package)
% 
% 
SEmex  = SolveEikonalmex(myTM,myBdy,myParam,myMetric);
tic
Mobj.distOB  = SEmex.Compute_Soln;
cd(CD)
if ftbverbose
    fprintf('end   : %s \n', subname)
end
return
