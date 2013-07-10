function ww3_to_swan_bndry();
	
% Generate SWAN boundary forcing by interpolating from WW3 output
%
% function swan2netcdf(matfile,ncfile,basename,first_time,last_time,increment);
%
% DESCRIPTION:
%    interpolate Hs,Tp,Dir from WW3 to a boundary forcing file for 
%    unstructured swan
%    this is an example file and will need to be modified for specific cases
%    Note that for the unstructured SWAN you can specify separate TPAR files
%    containing time series of hs,tp,dir for each node.  Also note that the nodes
%    are not necessarily nodes of the SWAN mesh.  They are specified in arclength
%    of the grid units from the first open boundary node (arclength 0).  SWAN
%    assembles a boundary segment by piecing together the nodes marked as boundary 
%    nodes in the node file (mark = 2).  This assumes somehow that the nodes are ordered
%    sequentially along the boundary arc which is in fact a major assumption.
%
%    Sample OBC section of a swan input file for the GoM domain is as follows where
%    15 points are used to specify the boundary forcing while the domain in fact has 
%    60 boundary points.  SWAN interpolates as necessary to force all the boundary nodes.
%    The large numbers are the arclengths in meters
%
% BOUNDSPEC SIDE 2 CLOCKWISE VARIABLE FILE &
%         0.00 'obc1.bnd' 1 &
%     52704.26 'obc2.bnd' 1 &
%    131926.06 'obc3.bnd' 1 &
%    255117.10 'obc4.bnd' 1 &
%    390381.71 'obc5.bnd' 1 &
%    559989.50 'obc6.bnd' 1 &
%    740759.98 'obc7.bnd' 1 &
%    924330.66 'obc8.bnd' 1 &
%   1104489.93 'obc9.bnd' 1 &
%   1295381.43 'obc10.bnd' 1 &
%   1480466.74 'obc11.bnd' 1 &
%   1641071.70 'obc12.bnd' 1 &
%   1750424.20 'obc13.bnd' 1 &
%   1828825.67 'obc14.bnd' 1 &
%   1951072.38 'obc15.bnd' 1

% INPUT 
%
% OUTPUT:
%   SWAN open boundary TPAR files obcXX.bnd
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================


swan_node_file = '../gom1/gom1.node';

% set year range
ybeg = 2007;
yend = 2007;

% increment for dumping forcing: 
% swan does not force at nodes, just at locations along the arclength
% of the boundary segment.  
inc = 4;

% matlab output file
outfile = ['swan_gom0_obc_2007' num2str(ybeg)];

% datadir - location of all ww3 grib files
ddir = '/Volumes/Data/Users/gcowles/Data/SCALLOP_RSA_2009/wavewatch_wna_grib/'; 

% read the swan node file and grab boundary nodes
[num,x,y,mark] = textread(swan_node_file,'%d %f %f %d\n','headerlines',1);
obc_nodes = find(mark==2);
nobc = prod(size(obc_nodes));
xtmp = x(obc_nodes);
ytmp = y(obc_nodes);
arc = zeros(nobc,1); 
for i=2:nobc
  arc(i) = arc(i-1) + sqrt( (xtmp(i)-xtmp(i-1))^2 + (ytmp(i)-ytmp(i-1))^2); 
end;  

% shift to discrete locations
pts = 1:inc:nobc;
pts(end) = nobc;
xobc = xtmp(pts); 
yobc = ytmp(pts); 
aobc = arc(pts); 
ndisc = prod(size(xobc));


% inverse project discrete locations to lon/lat 
%junk = 0;
fid = fopen('in.dat','w');
for i=1:ndisc  
  fprintf(fid,'%f %f\n',xobc(i),yobc(i));
end;
fclose(fid);
system('./project_cmd');
fid = fopen('out.dat','r');
for i=1:ndisc  
  C = textscan(fid, '%f %f', 1);
  lon_obc(i) = C{1};
  lat_obc(i) = C{2};
end;
fclose(fid);
fprintf('finished: projecting to lon/lat\n')


%---------------------------------------------------------
% read a sample grib file and reconstruct WW3 grid => xg,yg
%---------------------------------------------------------
fname = [ddir 'wna.hs.200711.grb'];
grib_struct=read_grib(fname,[1],'ScreenDiag',0);
xmin = grib_struct.gds.Lo1;
xmax = grib_struct.gds.Lo2;
ymin = grib_struct.gds.La1;
ymax = grib_struct.gds.La2;
il = grib_struct.gds.Ni;
jl = grib_struct.gds.Nj;
dlon = (xmax-xmin)/(il-1);
dlat = (ymax-ymin)/(jl-1);
lon = xmin:dlon:xmax;
lon = lon-360;
lat = ymin:dlat:ymax;  
[xg,yg] = meshgrid(lon,lat);


%---------------------------------------------------------
% extract Hs,Tp,Dir from each WWIII file 
%---------------------------------------------------------
ndays      = [31,28,31,30,31,30,31,31,30,31,30,31];
ndays_leap = [31,29,31,30,31,30,31,31,30,31,30,31];

% loop over data and count number of days of data to preallocate
icnt = 0;
for y=ybeg:yend  
for m=1:1 %debug 12;
  year = int2str(y);
  mnth = int2str(m);
  if(m>9)
    fname = [ddir 'wna.hs.' year mnth '.grb'];  
  else
    fname = [ddir 'wna.hs.' year '0' mnth '.grb'];  
  end;
  if(exist(fname)) 
    fprintf('file %s exists\n',fname)
    if(mod(y,4)==0) %leap year (note 2OOO was a leap year, 1900, 2100 are not)
      icnt = icnt + 8*ndays_leap(m);
    else
      icnt = icnt + 8*ndays(m);
    end; 
  else 
    fprintf('file %s does not exist\n',fname)
  end;
end;
end;
fprintf('number of frames in year %d\n',icnt);

% preallocate arrays 
hs = zeros(ndisc,icnt);
tp = zeros(ndisc,icnt);
dir = zeros(ndisc,icnt);
time = zeros(icnt,1);

%
hour = [0,3,6,9,12,15,18,21];

% read data into the arrays
icnt = 0;
for y=ybeg:yend;
for m=1:1 %debug 12;
  year = int2str(y);
  mnth = int2str(m);
  if(m>9)
    hsname = [ddir 'wna.hs.' year mnth '.grb'];
    tpname = [ddir 'wna.tp.' year mnth '.grb'];
    drname = [ddir 'wna.dp.' year mnth '.grb'];
  else
    hsname = [ddir 'wna.hs.' year '0' mnth '.grb'];
    tpname = [ddir 'wna.tp.' year '0' mnth '.grb'];
    drname = [ddir 'wna.dp.' year '0' mnth '.grb'];
  end;
  if(exist(hsname)) 
    fprintf('processing year %d month %d\n',y,m);
    if(mod(y,4)==0) %leap year
       
       nd = ndays_leap(m);
    else
       nd = ndays(m);
    end;
  
    n = 0;
    for d=1:nd  %day loop
    for l=1:8   %hour loop

      n = n + 1;
      icnt = icnt + 1;
  
      % read hs and interpolate onto boundary points 
      grib_struct=read_grib(hsname,[n],'ScreenDiag',0);
      var = reshape(grib_struct.fltarray,il,jl);
      hs(:,icnt) = interp2(xg,yg,var',lon_obc,lat_obc,'linear');
  
      %read tp from grib and interpolate onto boundary points 
      grib_struct=read_grib(tpname,[n],'ScreenDiag',0);
      var = reshape(grib_struct.fltarray,il,jl);
      tp(:,icnt) = interp2(xg,yg,var',lon_obc,lat_obc,'linear');

      %read dir from grib and interpolate onto boundary points 
      grib_struct=read_grib(drname,[n],'ScreenDiag',0);
      var = reshape(grib_struct.fltarray,il,jl);
      dir(:,icnt) = interp2(xg,yg,var',lon_obc,lat_obc,'linear');

      fprintf('processing time %s\n',grib_struct.stime);
  
      time(icnt) = greg2julian(y,m,d,hour(l),0,0);
      if(icnt>1);
      fprintf('processing frame %d %d %d %d %d %f\n',n,y,m,d,hour(l),time(icnt)-time(icnt-1));
       end;
    end;
    end;
  end; %file exists
  
end;
end;

% process data to fix boundaries where GOM open boundary is considered land in WaveWatch
for i=1:icnt
  hs(1:2,i) = hs(3,i);
  hs(end-1:end,i) = hs(end-2,i);
  tp(1:2,i) = tp(3,i);
  tp(end-1:end,i) = tp(end-2,i);
  dir(1:2,i) = dir(3,i);
  dir(end-1:end,i) = dir(end-2,i);
end;

% find points with NaN type data and use data from nearest neighbor 
pts = 1:ndisc;
ney = pts + 1;
ney2 = pts - 1;
ney(end) = pts(end)-1;
ney2(1)   = pts(1)+1;
obcs = pts;

for j=1:15
for i=1:icnt
  pts = find(hs(:,i)>99); hs(pts,i) = hs(ney(pts),i);
  pts = find(hs(:,i)>99); hs(pts,i) = hs(ney2(pts),i);
  pts = find(tp(:,i)>99); tp(pts,i) = tp(ney(pts),i);
  pts = find(tp(:,i)>99); tp(pts,i) = tp(ney2(pts),i);
  pts = find(dir(:,i)>1000); dir(pts,i) = dir(ney(pts),i);
  pts = find(dir(:,i)>1000); dir(pts,i) = dir(ney2(pts),i);
end;  
end;  

% dump data to a matlab object and save 
info = 'ww3 data interpolated onto gom0 open boundary';
info2 = 'time in julian day';
save(outfile,'info','info2','time','hs','tp','dir','lon_obc','lat_obc');  

figure
subplot(3,1,1)
for i=1:ndisc 
plot(hs(i,:)); hold on;
end;
subplot(3,1,2)
for i=1:ndisc
plot(tp(i,:)); hold on;
end;
subplot(3,1,3)
for i=1:ndisc
plot(dir(i,:)); hold on;
end;

%dump to separate swan forcing files 
for i=1:ndisc
  fname = ['obc' num2str(i) '.bnd'];
  fid = fopen(fname,'w');
  fid = fprintf(fid,'TPAR\n');
  for j=1:icnt
    [year,month,day,hour,mint,sec] = julian2greg(time(j));
    if(day < 10)
      daystr = ['0' int2str(day)];
    else
      daystr = int2str(day);
    end;
    if(month < 10)
      monthstr = ['0' int2str(month)];
    else
      monthstr = int2str(month);
    end;
    if(hour < 10)
      hourstr = ['0' int2str(hour)];
    else
      hourstr = int2str(hour);
    end;
    date = [' ' int2str(year) monthstr daystr '.' hourstr '00'];
    fprintf(fid,'%s %f %f %f %f\n',date,hs(i,j),tp(i,j),dir(i,j),10.);
  end;
  fclose(fid);
end;

% dump the main swan control file list
fname = 'cntrllist.txt';
fid = fopen(fname,'w');
for i=1:ndisc
  fname = ['"' 'obc' num2str(i) '.bnd' '"'];
  fprintf(fid,'%12.2f %s %d %s\n',aobc(i),fname,1,'&');
end;

