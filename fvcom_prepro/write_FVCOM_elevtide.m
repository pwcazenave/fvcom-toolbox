function write_FVCOM_elevtide(ObcNodes,JulianTime,SurfaceElevation,ElevationFile,MyTitle)
	
% Write an FVCOM surface elevation time series forcing file 
%
% function write_FVCOM_elevtide(ObcNodes,JulianTime,SurfaceElevation,SpectralFile,MyTitle)
%
% DESCRIPTION:
%    Write an FVCOM NetCDF surface elevation forcing file
%
% INPUT:
%   ObcNodes         = list of open boundary nodes of size [nObcs]
%   JulianTime       = list of modified Julian Dates of size [times] (but
%   defined as unlimited in the NetCDF file.
%   SurfaceElevation = list of surface elevation values of size [nObcs,
%   times]
%   ElevationFile    = name of NetCDF file
%   MyTitle          = case title, written as global attribute of NetCDF
%   file
%
% OUTPUT:
%    ElevationFile, A NetCDF FVCOM surface elevations tide forcing file
%
% EXAMPLE USAGE
%    write_FVCOM_elevtide(ObcNodes,JulianTime,SurfaceElevation,SpectralFile,MyTitle)
%
% Author(s):  
%    Pierre Cazenave (Plymouth Marine Laboratory)
% 
% Revision history
%    2012-08-08 First version.
%   
%==============================================================================

global ftbverbose 
report = false;
if(ftbverbose); report = true; end;
subname = 'write_FVCOM_elevtide';
if(report);  fprintf('\n'); end;
if(report); fprintf(['begin : ' subname '\n']); end;
%------------------------------------------------------------------------------
% Sanity check on input and dimensions
%------------------------------------------------------------------------------
nTimes = numel(JulianTime);
if(report); fprintf('Number of time steps %d\n',nTimes); end;

nObcs = numel(ObcNodes);
if(report); fprintf('Number of Open Boundary Nodes %d\n',nObcs); end;

[chk1,chk2] = size(SurfaceElevation);
if nObcs ~= chk1 || nTimes ~= chk2
    fprintf('Surface elevation dimensions do not match time series and number of boundary nodes.\n')
    fprintf('Surface elevation nodes and time sizes: (%d, %d)\n', chk1, chk2)
    fprintf('Boundary nodes size: %d\n', nObcs)
    fprintf('Times size: %d\n', nTimes)
	error('Input data sizes do not match. Check and try again.');
end

%%
%------------------------------------------------------------------------------
% Dump the file
%------------------------------------------------------------------------------

nc=netcdf.create(ElevationFile,'clobber');

% define global attributes
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type','FVCOM TIME SERIES ELEVATION FORCING FILE')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title',MyTitle)
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history','FILE CREATED using write_FVCOM_elevtide')

% define dimensions
nobc_dimid=netcdf.defDim(nc,'nobc',nObcs);
time_dimid=netcdf.defDim(nc,'time',netcdf.getConstant('NC_UNLIMITED'));
date_str_len_dimid=netcdf.defDim(nc,'DateStrLen',[26,1]);

% define variables and attributes
nobc_varid=netcdf.defVar(nc,'obc_nodes','NC_INT',nobc_dimid);
netcdf.putAtt(nc,nobc_varid,'long_name','Open Boundary Node Number');
netcdf.putAtt(nc,nobc_varid,'grid','obc_grid');

iint_varid=netcdf.defVar(nc,'iint','NC_INT',time_dimid);
netcdf.putAtt(nc,iint_varid,'long_name','internal mode iteration number');

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

Times_varid=netcdf.defVar(nc,'Times','NC_CHAR',[date_str_len_dimid, time_dimid]);
netcdf.putAtt(nc,Times_varid,'time_zone','UTC');

elevation_varid=netcdf.defVar(nc,'elevation','NC_FLOAT',[nobc_dimid, time_dimid]);
netcdf.putAtt(nc,elevation_varid,'long_name','Open Boundary Elevation');
netcdf.putAtt(nc,elevation_varid,'units','meters');

% end definitions
netcdf.endDef(nc);

% write data
netcdf.putVar(nc,nobc_varid,ObcNodes);
netcdf.putVar(nc,iint_varid,0,nTimes,1:nTimes);
netcdf.putVar(nc,time_varid,0,nTimes,JulianTime - 678942);
netcdf.putVar(nc,itime_varid,floor(JulianTime - 678942));
netcdf.putVar(nc,itime2_varid,0,nTimes,mod(JulianTime - 678942,1)*24*3600*1000);
nStringOut = char();
for i=1:nTimes
    nStringOut = [nStringOut, sprintf('%04i/%02i/%02i %02i:%02i:%02i       ',datevec(JulianTime(i)))];
end
netcdf.putVar(nc,Times_varid,nStringOut);
netcdf.putVar(nc,elevation_varid,SurfaceElevation);

% close file
netcdf.close(nc);

if(report); fprintf(['end   : ' subname '\n']); end;

