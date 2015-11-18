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
%   nverts = number of vertices in the mesh
%   tsl = temperature at level k (C)
%   ssl = salinity at level k (PSU)
%   filename  = filename to which to dump data
%   mytitle   = global attribute
%
% OUTPUT:
%    NetCDF file: filename
%
%
% EXAMPLE USAGE
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Chang Liu (University of Massachusetts Dartmouth)
%
% Revision history
%    2012-06-15 Added support for native MATLAB NetCDF routines. Requires
%    MATLAB 2010a or higher.
%
%==============================================================================

global ftbverbose

subname = 'write_FVCOM_obc_TS'
if ftbverbose
    fprintf('\nbegin : %s \n', subname);
end

% check dimensions
ksl = numel(zsl);

% if(numel(tsl) ~= ksl)
%   error('dimensions of ssl do not match zsl')
% end;
% if(numel(ssl) ~= ksl)
%   error('dimensions of ssl do not match zsl')
% end;

%------------------------------------------------------------------------------
% Dump to S/T profile to NetCDF file
%------------------------------------------------------------------------------
fprintf(['Dumping to NetCDF file: ',filename,'\n']);
fprintf('Size of T/S array: %d\n',ksl);

nc = netcdf.create(filename,'clobber');

netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title',mytitle)
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history', sprintf('File created with %s from the MATLAB fvcom-toolbox', subname))

% define dimensions
ksl_dimid=netcdf.defDim(nc,'ksl',ksl);
node_dimid=netcdf.defDim(nc,'node',nverts);
time_dimid=netcdf.defDim(nc,'time',netcdf.getConstant('NC_UNLIMITED'));

% define variables and attributes
time_varid=netcdf.defVar(nc,'time','NC_FLOAT',time_dimid);
netcdf.putAtt(nc,time_varid,'long_name','time');
netcdf.putAtt(nc,time_varid,'units','days since 0.0');
netcdf.putAtt(nc,time_varid,'time_zone','none');

itime_varid=netcdf.defVar(nc,'Itime','NC_INT',time_dimid);
netcdf.putAtt(nc,itime_varid,'units','days since 0.0');
netcdf.putAtt(nc,itime_varid,'time_zone','none');

itime2_varid=netcdf.defVar(nc,'Itime2','NC_INT',time_dimid);
netcdf.putAtt(nc,itime2_varid,'units','msec since 00:00:00');
netcdf.putAtt(nc,itime2_varid,'time_zone','none');

zsl_varid=netcdf.defVar(nc,'zsl','NC_FLOAT',ksl_dimid);
netcdf.putAtt(nc,zsl_varid,'long_name','standard z levels positive up');
netcdf.putAtt(nc,zsl_varid,'units','m');

ssl_varid=netcdf.defVar(nc,'ssl','NC_FLOAT',[node_dimid,ksl_dimid,time_dimid]);
netcdf.putAtt(nc,ssl_varid,'long_name','observed_salinity_profile');
netcdf.putAtt(nc,ssl_varid,'units','PSU');

tsl_varid=netcdf.defVar(nc,'tsl','NC_FLOAT',[node_dimid,ksl_dimid,time_dimid]);
netcdf.putAtt(nc,tsl_varid,'long_name','observed_temperature_profile');
netcdf.putAtt(nc,tsl_varid,'units','C');

% end definitions
netcdf.endDef(nc);


% write vars
netcdf.putVar(nc,time_varid,0,numel(time),time);
netcdf.putVar(nc,itime_varid,floor(time));
netcdf.putVar(nc,itime2_varid,0,numel(time),(mod(time,1)*24*3600*1000));
netcdf.putVar(nc,zsl_varid,zsl);
netcdf.putVar(nc,tsl_varid,tsl);
netcdf.putVar(nc,ssl_varid,ssl);


% for i=1:numel(time)
%     for k=1:ksl
%         nc{'tsl'}(i,k,:) = tsl(k);
%     end;
% end;
%
% for i=1:numel(time)
%     for k=1:ksl
%         nc{'ssl'}(i,k,:) = ssl(k);
%     end;
% end;
%
% Close the NetCDF file(s)
netcdf.close(nc);

if ftbverbose
    fprintf('end   : %s \n',  subname)
end



