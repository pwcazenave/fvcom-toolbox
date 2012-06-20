function ideal_swan_obc(SwanNodeFile,time,hs,tp,dir,inc,DirectSpread);
	
% Generate idealized SWAN boundary forcing
%
% function  ideal_swan_obc(SwanNodeFile,time,hs,tp,dir,inc,DirectSpread);
%
% DESCRIPTION:
%    interpolate Hs,Tp,Dir from WW3 from time series to unstructured SWAN forcing file
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
%   SwanNodeFile:   Swan node file (e.g. tst.node)
%   time:           time stamp of time series in modified Julian day
%   hs:             time series for significant wave height 
%   tp:             time series for peak period
%   dir:            time series for wave direction
%   inc:            dump TPAR obc forcing file every inc boundary points 
%   DirectSpread    directional spreading of incoming waves in degrees
%
% OUTPUT:
%   SWAN open boundary TPAR files obcXX.bnd
%   cntrllist.txt:  list of open boundary arclengths and obc file names
%                   this can essentially be pasted into the swan input file
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================


% read the swan node file and grab boundary nodes
[num,x,y,mark] = textread(SwanNodeFile,'%d %f %f %d\n','headerlines',1);
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
nTPAR = prod(size(xobc));


%---------------------------------------------------------
% set TPAR files for each open boundary forcing point
%---------------------------------------------------------


% preallocate arrays 
nTimes = prod(size(time));

%dump to separate swan forcing files 
for i=1:nTPAR
  fname = ['obc' num2str(i) '.bnd'];
  fid = fopen(fname,'w');
  fid = fprintf(fid,'TPAR\n');
  for j=1:nTimes
    [year,month,day,hour,mint,sec] = mjulian2greg(time(j));
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
    fprintf(fid,'%s %f %f %f %f\n',date,hs(j),tp(j),dir(j),DirectSpread);
  end;
  fclose(fid);
end;

 % dump the main swan control file list
 fname = 'cntrllist.txt';
 fid = fopen(fname,'w');
 for i=1:nTPAR
   fname = [' ''obc' num2str(i) '.bnd'' '];
   fprintf(fid,'%12.2f %s %d %s\n',aobc(i),fname,1,'&');
 end;

