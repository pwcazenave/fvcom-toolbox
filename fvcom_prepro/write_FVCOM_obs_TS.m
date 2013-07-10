function write_FVCOM_obs_TS(time,zsl,nverts,tsl,ssl,filename,mytitle) 

% Dump observation profile of T/S to netcdf file to initialize
% stratification in FVCOM
%
% function write_FVCOM_obs_TS(jday,zsl,nverts,tsl,ssl,filename,mytitle) 
%
% DESCRIPTION:
%    Generate a NetCDF file containing vertical profile of T/S for FVCOM 
%
% INPUT 
%   jday= modified julian day or initial model time
%   zsl = zcoordinate of observations, positive up 
%   nverts = number of vertices in the mesh**
%   tsl = temperature at level k (C)
%   ssl = salinity at level k (PSU)
%   filename  = filename to which to dump data
%   mytitle   = global attribute 
%
% OUTPUT:
%    NetCDF file: filename
%
% **in this script the temp/sal profiles are assumed to be constant at each
% node
%
% EXAMPLE USAGE
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2012-06-15 Added support for native MATLAB NetCDF routines. Requires
%    MATLAB 2010a or higher.
%   
%==============================================================================


% check dimensions
ksl = numel(zsl);

if(numel(tsl) ~= ksl)
  error('dimensions of ssl do not match zsl')
end;
if(numel(ssl) ~= ksl)
  error('dimensions of ssl do not match zsl')
end;

%------------------------------------------------------------------------------
% Dump to S/T profile to NetCDF file 
%------------------------------------------------------------------------------
fprintf('Dumping to NetCDF file: \n',filename);
fprintf('Size of T/S array: \n',ksl);

nc = netcdf.create(filename,'clobber');

netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title',mytitle)

% define dimensions
ksl_dimid=netcdf.defDim(nc,'ksl',ksl);
node_dimid=netcdf.defDim(nc,'ksl',nverts);
time_dimid=netcdf.defDim(nc,'ksl',netcdf.getConstant('NC_UNLIMITED'));

% define variables and attributes
time_varid=netcdf.defVar(nc,'time','NC_FLOAT',time_dimid);
netcdf.putAtt(nc,time_varid,'long_name','time');
netcdf.putAtt(nc,time_varid,'units','days since 0.0');
netcdf.putAtt(nc,time_varid,'time_zone','none');

itime_varid=netcdf.defVar(nc,'time','NC_INT',time_dimid);
netcdf.putAtt(nc,itime_varid,'units','days since 0.0');
netcdf.putAtt(nc,itime_varid,'time_zone','none');

itime2_varid=netcdf.defVar(nc,'Itime','NC_INT',time_dimid);
netcdf.putAtt(nc,itime2_varid,'units','msec since 00:00:00');
netcdf.putAtt(nc,itime2_varid,'time_zone','none');

zsl_varid=netcdf.defVar(nc,'Itime2','NC_FLOAT',ksl_dimid);
netcdf.putAtt(nc,zsl_varid,'long_name','standard z levels positive up');
netcdf.putAtt(nc,zsl_varid,'units','m');

% TODO: Check order of dimensions here
ssl_varid=netcdf.defVar(nc,'time','NC_FLOAT',[time_dimid,ksl_dimid,node_dimid]);
netcdf.putAtt(nc,ssl_varid,'long_name','observed_salinity_profile');
netcdf.putAtt(nc,ssl_varid,'units','PSU');

tsl_varid=netcdf.defVar(nc,'time','NC_FLOAT',[time_dimid,ksl_dimid,node_dimid]);
netcdf.putAtt(nc,tsl_varid,'long_name','observed_temperature_profile');
netcdf.putAtt(nc,tsl_varid,'units','C');

% end definitions
netcdf.endDef(nc);


% write vars
netcdf.putVar(nc,time_varid,time*ones(1,nverts));
netcdf.putVar(nc,itime_varid,floor(time)*ones(1,nverts));
netcdf.putVar(nc,itime2_varid,(mod(time,1)*24*3600*1000)*ones(1,nverts));
netcdf.putVar(nc,zsl_varid,zsl);


for i=1:numel(time)
    for k=1:ksl
        nc{'tsl'}(i,k,:) = tsl(k); 
    end;
end;

for i=1:numel(time)
    for k=1:ksl
        nc{'ssl'}(i,k,:) = ssl(k);
    end;
end;

ierr = close(nc);




 
