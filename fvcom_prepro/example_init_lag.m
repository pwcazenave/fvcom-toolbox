clear all; close all;

% example script: 
% initialize the online Lagrangian tracking for Julian day (realtime) forcing
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
warning off
type = 2;

% set mesh/bathymetry files
meshfile = '../Model_Maker/coos01_grd.dat';
bathfile = '../Model_Maker/coos01_dep.dat'; 
lagfile  = 'test_lag.nc';

% start/end time for particles
tbeg = greg2mjulian(2008,01,01,0,0,0); 
tend = greg2mjulian(2008,02,01,0,0,0); 

% read in mesh and bathymetry
Mobj = read_fvcom_mesh(meshfile);
Mobj.h = read_fvcom_bath(bathfile); Mobj.have_bath = true;
Mobj = setup_metrics(Mobj);


if(type==1) %initialize at all elements

xc   = Mobj.xc;
yc   = Mobj.yc;
nLag = Mobj.nElems;

elseif(type==2) %initialize along a line of interest
nLag = 10;  
p1 = [1.188363e6,194497];
p2 = [1.188548e6,194996];
xp = p1(1):(p2(1)-p1(1))/(nLag-1):p2(1);
yp = p1(2):(p2(2)-p1(2))/(nLag-1):p2(2);

end;


% plot to check
plot_field(Mobj,Mobj.h,'title','domain','withextra',false,'showgrid',false); hold on;
plot(xp,yp,'ro');


% dump the initial particle position file
nc = netcdf(lagfile,'clobber');
nc.references = 'http://fvcom.smast.umassd.edu';
nc.source = 'lag_init.m';
nc.info = 'debugging ';





  
% dimensions
nc('nparticles') = nLag;

% particle  vars
nc{'x'} = ncfloat('nparticles');
nc{'x'}.long_name = 'particle x position';
nc{'x'}.units = 'm'; 

nc{'y'} = ncfloat('nparticles');
nc{'y'}.long_name = 'particle y position';
nc{'y'}.units = 'm'; 

nc{'z'} = ncfloat('nparticles');
nc{'z'}.long_name = 'particle z position';
nc{'z'}.units = 'm'; 

nc{'pathlength'} = ncfloat('nparticles');
nc{'pathlength'}.long_name = 'particle integrated path length'; 
nc{'pathlength'}.units = 'm'; 

nc{'tbeg'} = ncfloat('nparticles');
nc{'tbeg'}.long_name = 'particle release time';
nc{'tbeg'}.units = 'days since 1858-11-17 00:00:00';
nc{'tbeg'}.format = 'modified julian day (MJD)';
nc{'tbeg'}.time_zone = 'UTC';

nc{'tend'} = ncfloat('nparticles');
nc{'tend'}.long_name = 'particle freeze time';
nc{'tend'}.units = 'days since 1858-11-17 00:00:00';
nc{'tend'}.format = 'modified julian day (MJD)';
nc{'tend'}.time_zone = 'UTC';


nc{'group'} = ncint('nparticles');
nc{'group'}.long_name = 'particle group'; 
nc{'group'}.units = '-'; 

nc{'mark'} = ncint('nparticles');
nc{'mark'}.long_name = 'particle mark'; 
nc{'mark'}.units = '-'; 

% dump vars
nc{'x'}(:) = xp;
nc{'y'}(:) = yp;
nc{'z'}(:) = 0.0;
nc{'tbeg'}(:) = tbeg;
nc{'tend'}(:) = tend;
nc{'group'}(:) = 1;
nc{'mark'}(:) = 0;
nc{'pathlength'}(:) = 0.0;
close(nc); 
