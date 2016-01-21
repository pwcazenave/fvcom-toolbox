function Mobj=distance_along_BC(Mobj,BCnodes)
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
cd /users/modellers/rito/matlab/HJB_Solver_Package/HJB_Solver_Package
% 
% 
SEmex  = SolveEikonalmex(myTM,myBdy,myParam,myMetric);
tic
Mobj.distOB  = SEmex.Compute_Soln;
cd(CD)
end

