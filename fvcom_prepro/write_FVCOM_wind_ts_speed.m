function write_FVCOM_wind_ts_speed(Mobj, WindFile, time, u10, v10)

% Write out time-varying/spatially constant wind forcing as speed.
%
% function write_FVCOM_wind_ts_speed(Mobj, WindFile, time, u10, v10)
%
% DESCRIPTION:
%    Write a time-varying, spatially constant wind file
%
% INPUT
%    Mobj - MATLAB mesh object
%    WindFile - output NetCDF filename (including path)
%    time - time in MJD
%    u10 - vector x component of wind field 10m above the surface.
%    v10 - vector y component of wind field 10m above the surface.
%
% Note: the shape of u10 and v10 must match that of time since this
% currently only outputs temporally varying wind (not spatially varying).
%
% OUTPUT:
%    NetCDF WindFile
%
% EXAMPLE USAGE
%    time = 0:0.25:31;
%    write_FVCOM_wind_ts_speed(...
%       'casename_wnd.nc',...
%       time, ones(size(time)),...
%       ones(size(time))*0.25);
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2012-10-08 Converted example_FVCOM_wind_ts_speed.m to use built-in
%   MATLAB NetCDF functions for creating the output file, eliminating the
%   need for the third party NetCDF library. Also added three additional
%   arguments to the function call (time and u and v vectors). u and v
%   vectors vary in time and space.
%
%==============================================================================
warning off
subname = 'example_FVCOM_wind_ts_speed';
global ftbverbose;
if(ftbverbose);
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end

nElems = Mobj.nElems;
nNodes = Mobj.nVerts;

%------------------------------------------------------------------------------
% write output to time and spatially-varying FVCOM wind file
%------------------------------------------------------------------------------

nc=netcdf.create(WindFile,'clobber');

% define global attributes
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'references','http://fvcom.smast.umassd.edu')
% netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'source','single-point time-dependent surface forcing')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'source','fvcom grid (unstructured) surface forcing')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'institution','Plymouth Marine Laboratory')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history','Generated using the fvcom-toolbox')

% dimensions
nele_dimid=netcdf.defDim(nc,'nele',nElems);
nvert_dimid=netcdf.defDim(nc,'node',nNodes);
time_dimid=netcdf.defDim(nc,'time',netcdf.getConstant('NC_UNLIMITED'));

% time vars
time_varid=netcdf.defVar(nc,'time','NC_FLOAT',time_dimid);
netcdf.putAtt(nc,time_varid,'long_name','time');
netcdf.putAtt(nc,time_varid,'units','days since 1858-11-17 00:00:00');
netcdf.putAtt(nc,time_varid,'time_zone','none');
netcdf.putAtt(nc,time_varid,'format','modified julian day (MJD)');

itime_varid=netcdf.defVar(nc,'Itime','NC_INT',time_dimid);
netcdf.putAtt(nc,itime_varid,'units','days since 1858-11-17 00:00:00');
netcdf.putAtt(nc,itime_varid,'time_zone','none');
netcdf.putAtt(nc,itime_varid,'format','modified julian day (MJD)');

itime2_varid=netcdf.defVar(nc,'Itime2','NC_INT',time_dimid);
netcdf.putAtt(nc,itime2_varid,'units','msec since 00:00:00');
netcdf.putAtt(nc,itime2_varid,'time_zone','none');

% Space and time variables
u10_varid=netcdf.defVar(nc,'U10','NC_FLOAT',[nele_dimid,time_dimid]);
netcdf.putAtt(nc,u10_varid,'long_name','Eastward Wind Velocity');
netcdf.putAtt(nc,u10_varid,'standard_name','Wind Velocity');
netcdf.putAtt(nc,u10_varid,'units','m/s');
netcdf.putAtt(nc,u10_varid,'type','data');

v10_varid=netcdf.defVar(nc,'V10','NC_FLOAT',[nele_dimid,time_dimid]);
netcdf.putAtt(nc,v10_varid,'long_name','Northward Wind Velocity');
netcdf.putAtt(nc,v10_varid,'standard_name','Wind Velocity');
netcdf.putAtt(nc,v10_varid,'units','m/s');
netcdf.putAtt(nc,v10_varid,'type','data');


% end definitions
netcdf.endDef(nc);

% dump time
netcdf.putVar(nc,time_varid,0,numel(time),time);
netcdf.putVar(nc,itime_varid,floor(time));
netcdf.putVar(nc,itime2_varid,0,numel(time),mod(time,1)*24*3600*1000);
netcdf.putVar(nc,u10_varid,[0,0],[nElems,numel(time)],u10);
netcdf.putVar(nc,v10_varid,[0,0],[nElems,numel(time)],v10);

% close file
netcdf.close(nc);

if(ftbverbose);
    fprintf(['end   : ' subname '\n'])
end



