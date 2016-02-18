% Example script: 
%
% Initialize the online Lagrangian tracking for Julian day (realtime) forcing.
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2016-02-18 Updated to use the function I created (init_lag) which was in 
%   turn based on the original version of this script.
%   
%=============================================================================

clear
close all

type = 1;

% Set mesh/bathymetry files
meshfile = './samples/tst_grd.dat';
bathfile = './samples_dep.dat'; 
lagfile  = 'test_lag.nc';

% start/end time for particles
tbeg = greg2mjulian(2008,01,01,0,0,0); 
tend = greg2mjulian(2008,02,01,0,0,0); 

% read in mesh and bathymetry
Mobj = read_fvcom_mesh(meshfile);
Mobj.h = read_fvcom_bath(bathfile); Mobj.have_bath = true;
Mobj = setup_metrics(Mobj);

if type==1 % initialize at all elements
    xc   = Mobj.xc;
    yc   = Mobj.yc;
    nLag = Mobj.nElems;
elseif type==2 % initialize along a line of interest
    nLag = 10;
    p1 = [1.188363e6,194497];
    p2 = [1.188548e6,194996];
    xp = p1(1):(p2(1)-p1(1))/(nLag-1):p2(1);
    yp = p1(2):(p2(2)-p1(2))/(nLag-1):p2(2);
end

% plot to check
plot_field(Mobj,Mobj.h,'title','domain','withextra',false,'showgrid',false); hold on;
plot(xp,yp,'ro');

% dump the initial particle position file
init_lag(Mobj, [tbeg, tend], lagfile)
