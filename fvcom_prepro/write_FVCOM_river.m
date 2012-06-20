function write_FVCOM_river(RiverFile,RiverName,nRivnodes,time,flux,temp,salt,RiverInfo1,RiverInfo2)
% write FVCOM 3.x NetCDF river file
%
% function write_FVCOM_river(RiverFile,RiverName,nRivnodes,time,flux,temp,salt,RiverInfo1,RiverInfo2)
%
% DESCRIPTION:
%    Write river flux, temperature, and salinity to an FVCOM river file
%    Note that it is assumed that the NetCDF file contains data for only
%    one river, even if it is split among multiple nodes.  The flux will be
%    set at each node as flux/nRivnodes where nRivnodes is the number of River
%    nodes.  Salinity and Temperature will be set the same at each node
%
% INPUT
%    RiverFile:   FVCOM 3.x NetCDF river forcing file
%    RiverName:   Name of the actual River
%    nRivnodes:   # of River nodes
%    time     :   timestamp in modified Julian day 
%    flux     :   Total river flux of same dimensions as time in m^3/s
%    temp     :   temperature in C of same dimensions as time
%    salt     :   salinity in PSU of same dimensions as time
%    RiverInfo1 : global attribute of file
%    RiverInfo2 : additional global attribute of file
%   
% OUTPUT:
%    FVCOM RiverFile with flux,temp,salt
%
% EXAMPLE USAGE
%  write_FVCOM_river('tst_riv.nc','Penobscot',3,time,flux,salt,'Penobscot Flux','source: USGS')  
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
warning off;

global ftbverbose;
if(ftbverbose);
subname = 'write_FVCOM_river';
fprintf('\n')
fprintf(['begin : ' subname '\n'])
end;


if(ftbverbose); 
  fprintf('creating river NetCDF file %s for River %s\n',RiverFile,RiverName); 
end;


nTimes = prod(size(flux));
if(ftbverbose);
  fprintf('# of river nodes: %d\n',nRivnodes);
  fprintf('# of time frames: %d\n',nTimes);
end;

[year,month,day,hour,mint,sec] = mjulian2greg(time(1));
if(ftbverbose); fprintf('river begins at: %d %d %d\n',year,month,day); end;
[year,month,day,hour,mint,sec] = mjulian2greg(time(end));
if(ftbverbose); fprintf('river ends at:   %d %d %d\n',year,month,day); end;

% set the flux
if(ftbverbose); fprintf('dividing flux into %d points\n',nRivnodes); end;
river_flux = zeros(nTimes,nRivnodes);
for i=1:nTimes
  river_flux(i,1:nRivnodes) = flux(i)/real(nRivnodes);
end;

% set temperature and salt
for i=1:nTimes
	river_salt(i,1:nRivnodes) = salt(i);
	river_temp(i,1:nRivnodes) = temp(i);
end;

% set some kind of sediment
coarse_sand = 15*ones(nTimes,nRivnodes);
medium_sand = 45*ones(nTimes,nRivnodes);
fine_sand   = 30*ones(nTimes,nRivnodes);


%--------------------------------------------------------------
% dump to netcdf file
%--------------------------------------------------------------

% open boundary forcing
nc = netcdf(RiverFile, 'clobber');       

nc.type = 'FVCOM RIVER FORCING FILE' ;
nc.title = RiverInfo1;   
nc.info =  RiverInfo2; 
nc.history = 'FILE CREATED using write_river_file.m' ;

% dimensions
nc('rivers') = nRivnodes; 
nc('namelen') = 26; 
nc('time') = 0; 

% variables
nc{'river_names'} = ncchar('rivers', 'namelen');

nc{'time'} = ncfloat('time');
nc{'time'}.long_name = 'time';  
nc{'time'}.units     = 'days since 0.0';  
nc{'time'}.time_zone = 'none';  

nc{'Itime'} = ncint('time');
nc{'Itime'}.units     = 'days since 0.0';  
nc{'Itime'}.time_zone = 'none';  

nc{'Itime2'} = ncint('time');
nc{'Itime2'}.units     = 'msec since 00:00:00';
nc{'Itime2'}.time_zone = 'none';  

nc{'river_flux'} = ncfloat('time','rivers');
nc{'river_flux'}.long_name = 'river runoff volume flux'; 
nc{'river_flux'}.units     = 'm^3s^-1';  

nc{'river_temp'} = ncfloat('time','rivers');
nc{'river_temp'}.long_name = 'river runoff temperature'; 
nc{'river_temp'}.units     = 'Celsius';  

nc{'river_salt'} = ncfloat('time','rivers');
nc{'river_salt'}.long_name = 'river runoff salinity'; 
nc{'river_salt'}.units     = 'PSU';  

% river names (must be 26 character strings)
for i=1:nRivnodes
  fname = [RiverName int2str(i)];
  temp  = '                          ';
  temp(1:length(fname)) = fname;
  nc{'river_names'}(i,:)   = temp;
end;

% dump dynamic data
for i=1:nTimes
  nc{'time'}(i) = time(i);
  nc{'Itime'}(i) = floor(time(i));
  nc{'Itime2'}(i) = mod(time(i),1)*24*3600*1000.;
  nc{'river_flux'}(i,1:nRivnodes) = river_flux(i,1:nRivnodes); 
  nc{'river_temp'}(i,1:nRivnodes) = river_temp(i,1:nRivnodes); 
  nc{'river_salt'}(i,1:nRivnodes) = river_salt(i,1:nRivnodes); 
  nc{'coarse_sand'}(i,1:nRivnodes) = coarse_sand(i,1:nRivnodes); 
  nc{'medium_sand'}(i,1:nRivnodes) = medium_sand(i,1:nRivnodes); 
  nc{'fine_sand'}(i,1:nRivnodes) = fine_sand(i,1:nRivnodes); 
end;

nc = close(nc);    


if(ftbverbose);
  fprintf(['end   : ' subname '\n'])
end;

