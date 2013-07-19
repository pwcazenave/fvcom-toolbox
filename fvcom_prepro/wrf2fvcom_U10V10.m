%function wrf2fvcom_U10V10(wrf_file,fvcom_grid_file,fvcom_forcing_file,infos)
%------------------------------------------------------------------------
% interpolate wind heat fields from WRF onto the FVCOM mesh 
%------------------------------------------------------------------------
%warning off
wrf_file = 'Skagit_WRF_U10V10_2009.nc';
fvcom_grid_file =a 'skg4.3_grd.dat';
fvcom_forcing_file = 'skagit_2009_U10V10.nc';
infos = 'wrf from D. Ralston';
%------------------------------------------------------------------------
% open wrf file and read header 
%------------------------------------------------------------------------
if(~exist(wrf_file))
   error(['file: ' wrf_file ' does not exist']);
end;

% open wrf data and check for time range
nc = netcdf(wrf_file);
tmp = nc{'Times'}(:);
wrftime = greg2mjulian(str2num(tmp(:,1:4)),str2num(tmp(:,6:7)),...
           str2num(tmp(:,9:10)),str2num(tmp(:,12:13)),0,0);
fprintf('beg time of WRF data %s\n',tmp(1,:));
fprintf('end time of WRF data %s\n',tmp(end,:));
ntimes = prod(size(wrftime));
wrflat = nc{'XLAT'}(:,:);
wrflon = nc{'XLONG'}(:,:);
[nlat,nlon] = size(wrflat);

wrfx = zeros(nlat,nlon);
wrfy = zeros(nlat,nlon);

% project wrf grid to Euclidean for interpolation
for i=1:nlat
  [wrflon(i,:),wrflat(i,:),wrfx(i,:),wrfy(i,:)] = ...
      my_project(wrflon(i,:),wrflat(i,:),wrfx(i,:),wrfy(i,:),'forward');
end;

% load FVCOM mesh
Mobj = read_fvcom_mesh(fvcom_grid_file);
tri = Mobj.tri;
x   = Mobj.x;
y   = Mobj.y;
nVerts = Mobj.nVerts;
nElems = Mobj.nElems;
nv   = tri;  
fprintf('info for fvcom domain\n');
fprintf('number of nodes: %d\n',nVerts);
fprintf('number of elems: %d\n',nElems);

xc = zeros(nElems,1);
yc = zeros(nElems,1);
for i=1:nElems
  xc(i) = sum(x(tri(i,1:3)))/3.;
  yc(i) = sum(y(tri(i,1:3)))/3.;
end;

%---------------------------------------------------------------
% dump header for netcdf FVCOM forcing file
%---------------------------------------------------------------
nc = netcdf(fvcom_forcing_file, 'clobber');            
nc.type = 'FVCOM U10/V10 Forcing File' ;
nc.source = 'fvcom grid (unstructured) surface forcing';
nc.references = 'http://fvcom.smast.umassd.edu, http://codfish.smast.umassd.edu'; 
nc.institution = 'School for Marine Science and Technology' ;
nc.history = 'wrf_2_fvcom_U10V10.m';
nc.infos = infos;
  
% dimensions
nc('three') = 3;
nc('nele') = int8(nElems);
nc('node') = nVerts;
nc('time') = 0;

% time vars
nc{'time'} = ncfloat('time');
nc{'time'}.long_name = 'time';
nc{'time'}.units = 'days since 1858-11-17 00:00:00';
nc{'time'}.format = 'modified julian day (MJD)';
nc{'time'}.time_zone = 'UTC';
  
nc{'Itime'} = ncint('time');
nc{'Itime'}.units = 'days since 1858-11-17 00:00:00';
nc{'Itime'}.format = 'modified julian day (MJD)';
nc{'Itime'}.time_zone = 'UTC';

nc{'Itime2'} = ncint('time');
nc{'Itime2'}.units = 'msec since 00:00:00';
nc{'Itime2'}.time_zone = 'UTC';

nc{'x'} = ncint('node');
nc{'x'}.long_name = 'nodal x-coordinate';
nc{'x'}.units = 'm';
  
nc{'y'} = ncint('node');
nc{'y'}.long_name = 'nodal y-coordinate';
nc{'y'}.units = 'm';

nc{'nv'} = ncint('three','nele');
nc{'nv'}.long_name = 'nodes surrounding element';
nc{'nv'}(1:3,1:nElems) = tri';

nc{'U10'} = ncfloat('time','nele');
nc{'U10'}.long_name = 'Eastward 10-m Velocity';
nc{'U10'}.standard_name = 'Eastward Wind Speed';
nc{'U10'}.units = 'm/s';
nc{'U10'}.grid = 'fvcom_grid';
nc{'U10'}.type = 'data';

nc{'V10'} = ncfloat('time','nele');
nc{'V10'}.long_name = 'Northward 10-m Velocity';
nc{'V10'}.standard_name = 'Northtward Wind Speed';
nc{'V10'}.units = 'm/s';
nc{'V10'}.grid = 'fvcom_grid';
nc{'V10'}.type = 'data';

nc{'U10_node'} = ncfloat('time','node');
nc{'U10_node'}.long_name = 'Eastward 10-m Velocity';
nc{'U10_node'}.standard_name = 'Eastward Wind Speed';
nc{'U10_node'}.units = 'm/s';
nc{'U10_node'}.grid = 'fvcom_grid';
nc{'U10_node'}.type = 'data';

nc{'V10_node'} = ncfloat('time','node');
nc{'V10_node'}.long_name = 'Northward 10-m Velocity';
nc{'V10_node'}.standard_name = 'Northtward Wind Speed';
nc{'V10_node'}.units = 'm/s';
nc{'V10_node'}.grid = 'fvcom_grid';
nc{'V10_node'}.type = 'data';


% dump time
nc{'time'}(1:ntimes) = wrftime;
nc{'Itime'}(1:ntimes) = floor(wrftime); 
nc{'Itime2'}(1:ntimes) = mod(wrftime,1)*24*3600*1000.;
nc{'x'}(1:nVerts) = x;
nc{'y'}(1:nVerts) = y;

% read data from WRF grid, interpolate to FVCOM mesh
fvcom_u10   = zeros(nElems,1);
fvcom_v10   = zeros(nElems,1);
fvcom_u10_node   = zeros(nVerts,1);
fvcom_v10_node   = zeros(nVerts,1);

nc2 = netcdf(wrf_file);
icnt = 1;
for i=1:ntimes
  fprintf('interpolating frame %d of %d\n',i,ntimes);
  U10  = nc2{'U10'}(i,:,:);
  V10  = nc2{'V10'}(i,:,:);
  fvcom_u10_node  = griddata(wrfx,wrfy,U10,x,y);
  fvcom_v10_node  = griddata(wrfx,wrfy,V10,x,y);
  for j=1:nElems
     fvcom_u10(j) = sum(fvcom_u10_node(tri(j,1:3)))/3.; 
     fvcom_v10(j) = sum(fvcom_v10_node(tri(j,1:3)))/3.; 
  end;
  nc{'U10'}(icnt,1:nElems) = fvcom_u10;
  nc{'V10'}(icnt,1:nElems) = fvcom_v10;
  nc{'U10_node'}(icnt,1:nVerts) = fvcom_u10_node;
  nc{'V10_node'}(icnt,1:nVerts) = fvcom_v10_node;
  icnt = icnt + 1;
end;
fprintf('interpolation complete\n');

ierr = close(nc);
ierr = close(nc2);
