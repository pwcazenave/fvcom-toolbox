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

ncep_u10_file = '/tmp/irish_sea/raw_data/uwnd.sig995.2006.nc';
ncep_v10_file = '/tmp/irish_sea/raw_data/vwnd.sig995.2006.nc';
fvcom_forcing_file = '/tmp/irish_sea/input/configs/irish_sea_v9/irish_sea_v9_wnd_ncep.nc';
infos = 'NCEP from ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface/';

if exist(ncep_u10_file, 'file') ~= 2
   error(['file: ' ncep_u10_file ' does not exist']);
end
if exist(ncep_v10_file, 'file') ~= 2
    error(['file: ' ncep_v10_file ' does not exist']);
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
nc_u10 = netcdf.open(ncep_u10_file, 'NOWRITE');
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

[nceplon, nceplat] = meshgrid(nceplonvector, nceplatvector);

[nlat,nlon] = size(nceplat);

ncepx = zeros(nlat,nlon);
ncepy = zeros(nlat,nlon);

% project NCEP grid to cartesian grid for interpolation. Use Zone 30 here
% because it's just about the centre of the world -6 to 0 degrees
% longitude.
for i=1:nlat
  [ncepx(i,:),ncepy(i,:)] = cart_project(nceplon(i,:),nceplat(i,:));
end

% Get the relevant bits from the FVCOM mesh object
tri = Mobj.tri;
x   = Mobj.x;
y   = Mobj.y;
nVerts = Mobj.nVerts;
nElems = Mobj.nElems;
nv   = tri;
if(ftbverbose);
    fprintf('info for FVCOM domain\n');
    fprintf('number of nodes: %d\n',nVerts);
    fprintf('number of elems: %d\n',nElems);
end

xc = zeros(nElems,1);
yc = zeros(nElems,1);
for i=1:nElems
  xc(i) = mean(x(tri(i,1:3)));
  yc(i) = mean(y(tri(i,1:3)));
end

%---------------------------------------------------------------
% Dump header for netcdf FVCOM forcing file
%---------------------------------------------------------------
nc = netcdf.create(fvcom_forcing_file, 'clobber');

netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type','FVCOM U10/V10 Forcing File')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'source','fvcom grid (unstructured) surface forcing')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'references','http://fvcom.smast.umassd.edu, http://codfish.smast.umassd.edu')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'institution','Plymouth Marine Laboratory')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'institution','ncep_2_fvcom_U10V10.m')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'infos',infos)

% dimensions
three_dimid=netcdf.defDim(nc,'three',3);
nele_dimid=netcdf.defDim(nc,'nele',nElems);
node_dimid=netcdf.defDim(nc,'node',nVerts);
time_dimid=netcdf.defDim(nc,'time',netcdf.getConstant('NC_UNLIMITED'));

% time vars
time_varid=netcdf.defVar(nc,'time','NC_FLOAT',time_dimid);
netcdf.putAtt(nc,time_varid,'long_name','time');
netcdf.putAtt(nc,time_varid,'units','days since 1858-11-17 00:00:00');
netcdf.putAtt(nc,time_varid,'format','modified julian day (MJD)');
netcdf.putAtt(nc,time_varid,'time_zone','UTC');

itime_varid=netcdf.defVar(nc,'Itime','NC_INT',time_dimid);
netcdf.putAtt(nc,itime_varid,'units','days since 1858-11-17 00:00:00');
netcdf.putAtt(nc,itime_varid,'format','modified julian day (MJD)');
netcdf.putAtt(nc,itime_varid,'time_zone','UTC');

itime2_varid=netcdf.defVar(nc,'Itime2','NC_INT',time_dimid);
netcdf.putAtt(nc,itime2_varid,'units','msec since 00:00:00');
netcdf.putAtt(nc,itime2_varid,'time_zone','UTC');

% space vars
x_varid=netcdf.defVar(nc,'x','NC_FLOAT',node_dimid);
netcdf.putAtt(nc,x_varid,'long_name','nodal x-coordinate');
netcdf.putAtt(nc,x_varid,'units','m');

y_varid=netcdf.defVar(nc,'y','NC_FLOAT',node_dimid);
netcdf.putAtt(nc,y_varid,'long_name','nodal y-coordinate');
netcdf.putAtt(nc,y_varid,'units','m');

nv_varid=netcdf.defVar(nc,'nv','NC_FLOAT',[nele_dimid, three_dimid]);
netcdf.putAtt(nc,nv_varid,'long_name','nodes surrounding element');

% wind vars
u10_varid=netcdf.defVar(nc,'U10','NC_FLOAT',[nele_dimid, time_dimid]);
netcdf.putAtt(nc,u10_varid,'long_name','Eastward 10-m Velocity');
netcdf.putAtt(nc,u10_varid,'standard_name','Eastward Wind Speed');
netcdf.putAtt(nc,u10_varid,'units','m/s');
netcdf.putAtt(nc,u10_varid,'grid','fvcom_grid');
netcdf.putAtt(nc,u10_varid,'type','data');

v10_varid=netcdf.defVar(nc,'V10','NC_FLOAT',[nele_dimid, time_dimid]);
netcdf.putAtt(nc,v10_varid,'long_name','Northward 10-m Velocity');
netcdf.putAtt(nc,v10_varid,'standard_name','Northward Wind Speed');
netcdf.putAtt(nc,v10_varid,'units','m/s');
netcdf.putAtt(nc,v10_varid,'grid','fvcom_grid');
netcdf.putAtt(nc,v10_varid,'type','data');

u10_node_varid=netcdf.defVar(nc,'U10_node','NC_FLOAT',[node_dimid, time_dimid]);
netcdf.putAtt(nc,u10_node_varid,'long_name','Eastward 10-m Velocity');
netcdf.putAtt(nc,u10_node_varid,'standard_name','Eastward Wind Speed');
netcdf.putAtt(nc,u10_node_varid,'units','m/s');
netcdf.putAtt(nc,u10_node_varid,'grid','fvcom_grid');
netcdf.putAtt(nc,u10_node_varid,'type','data');

v10_node_varid=netcdf.defVar(nc,'V10_node','NC_FLOAT',[node_dimid, time_dimid]);
netcdf.putAtt(nc,v10_node_varid,'long_name','Northward 10-m Velocity');
netcdf.putAtt(nc,v10_node_varid,'standard_name','Northward Wind Speed');
netcdf.putAtt(nc,v10_node_varid,'units','m/s');
netcdf.putAtt(nc,v10_node_varid,'grid','fvcom_grid');
netcdf.putAtt(nc,v10_node_varid,'type','data');

% end definitions
netcdf.endDef(nc);

% write the triangulation data (nv)
netcdf.putVar(nc, nv_varid, tri');

% write the times
netcdf.putVar(nc,time_varid,0,ntimes,nceptime);
netcdf.putVar(nc,itime_varid,0,ntimes,floor(nceptime));
netcdf.putVar(nc,itime2_varid,0,ntimes,mod(nceptime,1)*24*3600*1000);
netcdf.putVar(nc,x_varid,x);
netcdf.putVar(nc,y_varid,y);


% read data from NCEP grid, interpolate to FVCOM mesh
fvcom_u10   = zeros(nElems,ntimes);
fvcom_v10   = zeros(nElems,ntimes);
fvcom_u10_node   = zeros(nVerts,ntimes);
fvcom_v10_node   = zeros(nVerts,ntimes);

% Open the as yet unopened V10 NCEP NetCDF file
nc_v10 = netcdf.open(ncep_v10_file, 'NOWRITE');

% Find the necessary variables
u10_varid_NCEP = netcdf.inqVarID(nc_u10, 'uwnd');
v10_varid_NCEP = netcdf.inqVarID(nc_v10, 'vwnd');

% The NCEP data are packed as integers. The following equation describes
% how to unpack them
%     unpacked value = add_offset + ( (packed value) * scale_factor )
% 
scale_factor = netcdf.getAtt(nc_u10,u10_varid_NCEP,'scale_factor','single');
add_offset = netcdf.getAtt(nc_u10,u10_varid_NCEP,'add_offset','single');

% Get the U10 and V10 data
U10 = netcdf.getVar(nc_u10, u10_varid_NCEP, 'single');
V10 = netcdf.getVar(nc_v10, v10_varid_NCEP, 'single');

% Scale by the factor
U10 = double(add_offset + (U10.*scale_factor));
V10 = double(add_offset + (V10.*scale_factor));

for i=1:ntimes
    fprintf('interpolating frame %d of %d\n', i, ntimes);

    fvcom_u10_node(:,i) = griddata(ncepx,ncepy,U10(:,:,i)',x,y);
    fvcom_v10_node(:,i) = griddata(ncepx,ncepy,V10(:,:,i)',x,y);
    for j=1:nElems
     fvcom_u10(j,i) = mean(fvcom_u10_node(tri(j,1:3))); 
     fvcom_v10(j,i) = mean(fvcom_v10_node(tri(j,1:3))); 
    end
%     nc{'U10'}(i,1:nElems) = fvcom_u10;
%     nc{'V10'}(i,1:nElems) = fvcom_v10;
%     nc{'U10_node'}(i,1:nVerts) = fvcom_u10_node;
%     nc{'V10_node'}(i,1:nVerts) = fvcom_v10_node;
end
fprintf('interpolation complete\n');

% Write the data to the FVCOM NetCDF input file
netcdf.putVar(nc, u10_varid, [0, i], [nElems, ntimes], fvcom_u10);
netcdf.putVar(nc, v10_varid, [0, i], [nElems, ntimes], fvcom_v10);
netcdf.putVar(nc, u10_node_varid, [0, i], [nVerts, ntimes], fvcom_u10_node);
netcdf.putVar(nc, v10_node_varid, [0, i], [nVerts, ntimes], fvcom_v10_node);

% Close the output file
netcdf.close(nc);
% Close the NCEP NetCDF files
netcdf.close(nc2);
netcdf.close(ncep_u10_file);
netcdf.close(ncep_v10_file);

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

function [out_east,out_north] = cart_project(in_east,in_north)
% Custom version of the my_project function which is provided with
% fvcom-toolbox. This one only goes from lat/long to eastings and
% northings. It also assumes zone 30N for the transformation since it's
% just about central (-6W to 0).
m_proj('UTM','longitude',[-180,-180],'latitude',[-90,90],'zone',30,'hemisphere','north','ellipsoid','wgs84')
[out_east, out_north] = m_ll2xy(in_east, in_north, 'clip', 'off');

