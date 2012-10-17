function [fvcom_u10_node, fvcom_v10_node] = ncep2fvcom_U10V10(Mobj,ncepx,ncepy,U10,V10,nceptime,fvcom_forcing_file,infos)
% Interpolate NCEP reanalysis wind speed data onto a given FVCOM grid
%
% ncep2fvcom_U10V10(Mobj,ncepx,ncepy,U10,V10,nceptime,fvcom_forcing_file,infos)
% 
% DESCRIPTION:
%   Takes a given NCEP reanalysis grid file and interpolates the U10 and
%   V10 values onto the specified FVCOM grid file. 
%   
% INPUT:
%   Mobj - MATLAB mesh object
%   ncepx - x data (probably best in cartesian for the interpolation)
%   ncepy - data (probably best in cartesian for the interpolation)
%   U10 - u-component wind data
%   V10 - v-component wind data
%   nceptime - NCEP Time vector (in Modified Julian Days)
%   fvcom_forcing_file - FVCOM forcing file name
%   infos - Additional remarks to be written to the "infos" NetCDF variable
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
%   2012-10-16 Removed the code to read the NCEP file. Instead, farmed that
%   out to a new function (read_NCEP_wind) so that the relevant section can
%   be more readily extracted (rather than using the entire globe's data:
%   it's easier to subsample and provide the subsampled data here). 
% 
%==========================================================================

warning off

if nargin ~= 8
    error('Incorrect number of arguments')
end

subname = 'ncep2fvcom_U10V10';

global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

%--------------------------------------------------------------------------
% Get the relevant bits from the FVCOM mesh object
%--------------------------------------------------------------------------
tri = Mobj.tri;
x   = Mobj.x;
y   = Mobj.y;
nVerts = Mobj.nVerts;
nElems = Mobj.nElems;
if(ftbverbose);
    fprintf('info for FVCOM domain\n');
    fprintf('number of nodes: %d\n',nVerts);
    fprintf('number of elems: %d\n',nElems);
end

xc = nodes2elems(x, Mobj);
yc = nodes2elems(y, Mobj);

ntimes = numel(nceptime);

% Interpolate NCEP data to FVCOM mesh
fvcom_u10   = zeros(nElems,ntimes);
fvcom_v10   = zeros(nElems,ntimes);
fvcom_u10_node   = zeros(nVerts,ntimes);
fvcom_v10_node   = zeros(nVerts,ntimes);

for i=1:ntimes
    fprintf('interpolating frame %d of %d\n', i, ntimes);

    fvcom_u10_node(:,i) = griddata(ncepx,ncepy,U10(:,:,i),x,y);
    fvcom_v10_node(:,i) = griddata(ncepx,ncepy,V10(:,:,i),x,y);
    for j=1:nElems
     fvcom_u10(j,i) = mean(fvcom_u10_node(tri(j,1:3))); 
     fvcom_v10(j,i) = mean(fvcom_v10_node(tri(j,1:3))); 
    end
end
fprintf('interpolation complete\n');


%--------------------------------------------------------------------------
% Dump header and data for netcdf FVCOM forcing file
%--------------------------------------------------------------------------
nc = netcdf.create(fvcom_forcing_file, 'clobber');

netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type','FVCOM U10/V10 Forcing File')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'source','fvcom grid (unstructured) surface forcing')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'references','http://fvcom.smast.umassd.edu, http://codfish.smast.umassd.edu')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'institution','Plymouth Marine Laboratory')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'institution','ncep_2_fvcom_U10V10.m')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'infos',infos)

% Dimensions
three_dimid=netcdf.defDim(nc,'three',3);
nele_dimid=netcdf.defDim(nc,'nele',nElems);
node_dimid=netcdf.defDim(nc,'node',nVerts);
time_dimid=netcdf.defDim(nc,'time',netcdf.getConstant('NC_UNLIMITED'));

% Time vars
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

% Space vars
x_varid=netcdf.defVar(nc,'x','NC_FLOAT',node_dimid);
netcdf.putAtt(nc,x_varid,'long_name','nodal x-coordinate');
netcdf.putAtt(nc,x_varid,'units','m');

y_varid=netcdf.defVar(nc,'y','NC_FLOAT',node_dimid);
netcdf.putAtt(nc,y_varid,'long_name','nodal y-coordinate');
netcdf.putAtt(nc,y_varid,'units','m');

nv_varid=netcdf.defVar(nc,'nv','NC_FLOAT',[nele_dimid, three_dimid]);
netcdf.putAtt(nc,nv_varid,'long_name','nodes surrounding element');

% Wind vars
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

% End definitions
netcdf.endDef(nc);

% Write data
netcdf.putVar(nc,nv_varid, tri');
netcdf.putVar(nc,time_varid,0,ntimes,nceptime);
netcdf.putVar(nc,itime_varid,0,ntimes,floor(nceptime));
netcdf.putVar(nc,itime2_varid,0,ntimes,mod(nceptime,1)*24*3600*1000);
netcdf.putVar(nc,x_varid,x);
netcdf.putVar(nc,y_varid,y);
netcdf.putVar(nc,u10_varid,[0,0],[nElems,ntimes],fvcom_u10);
netcdf.putVar(nc,v10_varid,[0,0],[nElems,ntimes],fvcom_v10);
netcdf.putVar(nc,u10_node_varid,[0,0],[nVerts,ntimes],fvcom_u10_node);
netcdf.putVar(nc,v10_node_varid,[0,0],[nVerts,ntimes],fvcom_v10_node);

% Close the NetCDF files
netcdf.close(nc);
netcdf.close(ncep_u10_file);
netcdf.close(ncep_v10_file);
