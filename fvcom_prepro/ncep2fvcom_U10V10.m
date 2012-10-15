function ncep2fvcom_U10V10(Mobj,ncep_u10_file,ncep_v10_file,fvcom_forcing_file,infos)
% Interpolate NCEP reanalysis wind speed data onto a given FVCOM grid
%
% function ncep2fvcom_U10V10()
% 
% DESCRIPTION:
%   Takes a given NCEP reanalysis grid file and interpolates the U10 and
%   V10 values onto the specified FVCOM grid file. 
%   
% INPUT:
%   Mobj - MATLAB mesh object
%   NCEP NetCDF U10 filename (and path)
%   NCEP NetCDF V10 filename (and path)
%   FVCOM forcing file name
%   Additional remarks to be written to the "infos" NetCDF variable
% 
% OUTPUT:
%   FVCOM wind speed forcing file
% 
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
% 
% Revision history:
%   2012-10-15 First version based on ncep2fvcom_U10V10.m in the
%   fvcom-toolbox.
% 
%==========================================================================

warning off

if nargin ~= 5
    error('Incorrect number of arguments')
end

subname = 'ncep2fvcom_U10V10';

global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

ncep_file = '/tmp/irish_sea/raw_data/uwnd.sig995.2006.nc';
fvcom_grid_file = '/data/medusa/pica/models/FVCOM/irish_sea/input/configs/irish_sea_v9/irish_sea_v9_grd.dat';
fvcom_forcing_file = '/tmp/irish_sea/input/configs/irish_sea_v9/irish_sea_v9_wnd_ncep.nc';
infos = 'NCEP from ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface/';

if exist(ncep_file, 'file') ~= 2
   error(['file: ' ncep_file ' does not exist']);
end

%--------------------------------------------------------------------------
% Open NCEP data and check for time range
%--------------------------------------------------------------------------

% Get the year from the NCEP file name
ncep_u10_year = getYear(ncep_u10_file);
ncep_v10_year = getYear(ncep_v10_file);
% Check both files are from the same year
if (ncep_u10_year - ncep_v10_year) ~= 0
    error('Input U and V wind data files are from different years')
end

% Get the time. It's stored relative to the 1st January for the given year.
% We're assuming that since both files are for the same year, we don't have
% to pull the 'time' variable out from both (they should be identical).
% Currently uses the U vector.
nc_u10 = netcdf.open(ncep_file, 'NOWRITE');
time_varid = netcdf.inqVarID(nc_u10, 'time');
nceptimehours = netcdf.getVar(nc_u10, time_varid);

% NCEP dates are relative to 0001/01/01 00:00:00 and stored in hours.
% MATLAB's dates are relative to 0000/00/00 00:00:00 and stored in days.
% Need to add a year and a day to the NCEP time when converting.
nceptimedays = datevec(nceptimehours/24 + (datenum(1, 0, 0) - 1));
nceptime = greg2mjulian(nceptimedays(:,1), nceptimedays(:,2),...
    nceptimedays(:,3), nceptimedays(:,4), nceptimedays(:,5),...
    nceptimedays(:,1));

if(ftbverbose);
    fprintf('beg time of NCEP data %04i/%02i/%02i %02i:%02i:%02i\n',nceptimedays(1,:));
    fprintf('end time of NCEP data %04i/%02i/%02i %02i:%02i:%02i\n',nceptimedays(end,:));
end

ntimes = numel(nceptime);

% Get the geographical information from the NCEP data. Again, use the U10
% file only (we're assuming they're both global).
lat_varid = netcdf.inqVarID(nc_u10, 'lat');
lon_varid = netcdf.inqVarID(nc_u10, 'lon');
nceplatvector = netcdf.getVar(nc_u10, lat_varid);
nceplonvector = netcdf.getVar(nc_u10, lon_varid);

[nceplat, nceplon] = meshgrid(nceplonvector, nceplatvector);

[nlat,nlon] = size(nceplat);

ncepx = zeros(nlat,nlon);
ncepy = zeros(nlat,nlon);

% project NCEP grid to cartesian grid for interpolation. Use Zone 30 here
% because it's just about the centre of the world -6 to 0 degrees
% longitude.
utmZones=cellfun(@(x) repmat(x,length(Mobj.x),1),'30 U','uni',false);
for i=1:nlat
  [nceplon(i,:),nceplat(i,:),ncepx(i,:),ncepy(i,:)] = ...
      my_project(nceplon(i,:),nceplat(i,:),ncepx(i,:),ncepy(i,:),'forward');
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
nc.history = 'ncep_2_fvcom_U10V10.m';
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
nc{'time'}(1:ntimes) = nceptime;
nc{'Itime'}(1:ntimes) = floor(nceptime); 
nc{'Itime2'}(1:ntimes) = mod(nceptime,1)*24*3600*1000.;
nc{'x'}(1:nVerts) = x;
nc{'y'}(1:nVerts) = y;

% read data from NCEP grid, interpolate to FVCOM mesh
fvcom_u10   = zeros(nElems,1);
fvcom_v10   = zeros(nElems,1);
fvcom_u10_node   = zeros(nVerts,1);
fvcom_v10_node   = zeros(nVerts,1);

nc2 = netcdf(ncep_file);
icnt = 1;
for i=1:ntimes
  fprintf('interpolating frame %d of %d\n',i,ntimes);
  U10  = nc2{'U10'}(i,:,:);
  V10  = nc2{'V10'}(i,:,:);
  fvcom_u10_node  = griddata(ncepx,ncepy,U10,x,y);
  fvcom_v10_node  = griddata(ncepx,ncepy,V10,x,y);
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

function year = getYear(file)
% Extract the year from a give NCEP file name (either 'uwnd.sig995.YYYY.nc'
% or 'vwnd.sig995.YYYY.nc'). Files are those downloaded from
% ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface/
% 
[~, tmp_year, ~] = fileparts(file);
tmp_year = regexp(tmp_year, '\.', 'split');
year = str2double(tmp_year(end));
if ~isnumeric(tmp_year)
    error('Could not parse the NCEP year from the NCEP file name. Expecting ''uwnd.sig995.YYYY.nc'' and ''vwnd.sig995.YYYY.nc''')
end
