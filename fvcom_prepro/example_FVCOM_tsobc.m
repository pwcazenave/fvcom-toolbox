function example_FVCOM_tsobc(basename,time,nSiglay)
% example file for dumping a file to force temperature and salinity at the open b.
%
% function example_FVCOM_tsobc()
%
% DESCRIPTION:
%    Setup a sample FVCOM hydrographic open boundary forcing file
%
% INPUT
%    Model case name
%    Time
%    Number of sigma layers
%    Number of sigma levels
%
% OUTPUT:
%    FVCOM hydrographic open boundary file
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
% 
% Revision history
%    2012-06-15 Added support for native MATLAB NetCDF routines. Requires
%    MATLAB 2010a or higher.
%    2012-07-16 Removed hard-coded nSiglay and nSiglev and instead moved to
%    arguments list.
% 
%==============================================================================

warning off;


subname = 'example_FVCOM_tsobc';
global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

fvcom_bathy = [basename, '_dep.dat'];
fvcom_obc   = [basename, '_obc.dat'];
tsOBCFile = [basename, '_tsobc.nc'];

%------------------------------------------------------------------------------
% read in the FVCOM open boundary node data (need node numbers and dimension)
%------------------------------------------------------------------------------
fid = fopen(fvcom_obc,'r');
if(fid  < 0)
  error(['file: ' fvcom_obc ' does not exist']);
end;
C = textscan(fid, '%s %s %s %s %d', 1);
nObc = C{5};
obc_nodes = zeros(nObc,1);
if(ftbverbose); fprintf('reading obc file\n'); end;
if(ftbverbose); fprintf('# nodes %d\n',nObc); end;
for i=1:nObc
  C = textscan(fid, '%d %d %d', 1);
  obc_nodes(i) = C{2};
end;

if(ftbverbose); fprintf('obc reading complete\n');end;

%------------------------------------------------------------------------------
% read in the FVCOM bathymetry data (need bathymetry on open boundary nodes)
%------------------------------------------------------------------------------
fid = fopen(fvcom_bathy,'r');
if(fid  < 0)
  error(['file: ' fvcom_bathy ' does not exist']);
end;
C = textscan(fid, '%s %s %s %d', 1);
Nverts = C{4};
h = zeros(Nverts,1);
if(ftbverbose); fprintf('reading bathymetry file\n');end;
if(ftbverbose); fprintf('# nodes %d\n',Nverts);end;
for i=1:Nverts
  C = textscan(fid, '%f %f %f', 1);
  h(i) = C{3};
end;
if(ftbverbose); fprintf('min depth %f max depth %f\n',min(h),max(h));end;
if(ftbverbose); fprintf('bathymetry reading complete\n');end;
fclose(fid);

%--------------------------------------------------------------
% set variables for NetCDF file
%--------------------------------------------------------------

% extract bathymetry at open boundary nodes
obc_h = h(obc_nodes);

% time
% time = 0:1:31.;
nTimes = numel(time);

% set siglev/siglay
% nSiglay = 10;
% nSiglev = 11;
nSiglev = nSiglay + 1;
inc = 1./real(nSiglay);
siglev = 0:-inc:-1;
for i=1:nSiglay
	siglay(i) = mean(siglev(i:i+1));
end;


% initialize temperature/salinity arrays
temp = zeros(nObc,nSiglay,nTimes);
salt = zeros(nObc,nSiglay,nTimes);

% set variable temperature and salinity
% for i=1:nTimes
% 	obc_temp(i) = 18. + 2.*real(i-1)/nTimes;
% 	obc_salt(i) = 30. - 5.*real(i-1)/nTimes;
% end

% set a constant temperature and salinity
obc_temp = ones(1,nTimes)*13;
obc_salt = ones(1,nTimes)*35;

%--------------------------------------------------------------
% dump to netcdf file
%--------------------------------------------------------------

% open boundary forcing
nc = netcdf.create(tsOBCFile, 'clobber');

% define global attributes
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type','FVCOM RIVER FORCING FILE')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title','simple open boundary hydrography test')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type','FVCOM TIME SERIES OBC TS FILE')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history','generated using the fvcom-toolbox')


% define dimensions
nobc_dimid=netcdf.defDim(nc,'nobc',nObc);
datestrlen_dimid=netcdf.defDim(nc,'Datestrln',26);
time_dimid=netcdf.defDim(nc,'time',netcdf.getConstant('NC_UNLIMITED'));
siglay_dimid=netcdf.defDim(nc,'siglay',nSiglay);
siglev_dimid=netcdf.defDim(nc,'siglev',nSiglev);

% variables
% nc{'river_names'} = ncchar('rivers', 'namelen');

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

nobc_varid=netcdf.defVar(nc,'obc_nodes','NC_INT',nobc_dimid);
netcdf.putAtt(nc,nobc_varid,'long_name','Open Boundary Node Number');
netcdf.putAtt(nc,nobc_varid,'grid','obc_grid');

obc_h_varid=netcdf.defVar(nc,'obc_h','NC_FLOAT',nobc_dimid);
netcdf.putAtt(nc,obc_h_varid,'long_name','ocean boundary depth');
netcdf.putAtt(nc,obc_h_varid,'units','m');
netcdf.putAtt(nc,obc_h_varid,'grid','obc_grid');

obc_siglev_varid=netcdf.defVar(nc,'siglev','NC_FLOAT',siglev_dimid);
netcdf.putAtt(nc,obc_siglev_varid,'long_name','ocean_sigma/general_coordinate');
netcdf.putAtt(nc,obc_siglev_varid,'grid','obc_grid');

obc_siglay_varid=netcdf.defVar(nc,'siglay','NC_FLOAT',siglay_dimid);
netcdf.putAtt(nc,obc_siglay_varid,'long_name','ocean_sigma/general_coordinate');
netcdf.putAtt(nc,obc_siglay_varid,'grid','obc_grid');

obc_temp_varid=netcdf.defVar(nc,'obc_temp','NC_FLOAT',[nobc_dimid,siglay_dimid,time_dimid]);
netcdf.putAtt(nc,obc_temp_varid,'long_name','sea_water_temperature');
netcdf.putAtt(nc,obc_temp_varid,'units','Celcius');
netcdf.putAtt(nc,obc_temp_varid,'grid','obc_grid');

obc_salinity_varid=netcdf.defVar(nc,'obc_salinity','NC_FLOAT',[nobc_dimid,siglay_dimid,time_dimid]);
netcdf.putAtt(nc,obc_salinity_varid,'long_name','sea_water_salinity');
netcdf.putAtt(nc,obc_salinity_varid,'units','PSU');
netcdf.putAtt(nc,obc_salinity_varid,'grid','obc_grid');

% end definitions
netcdf.endDef(nc);

% write data
netcdf.putVar(nc,nobc_varid,obc_nodes);
netcdf.putVar(nc,obc_h_varid,obc_h);
netcdf.putVar(nc,obc_siglev_varid,siglev);
netcdf.putVar(nc,obc_siglay_varid,siglay);
netcdf.putVar(nc,time_varid,0,numel(time),time);
netcdf.putVar(nc,itime_varid,floor(time));
netcdf.putVar(nc,itime2_varid,0,numel(time),mod(time,1)*24*3600*1000);

% Create 3D array from three 1D arrays
for i=1:nObc
    for j=1:nSiglay
        temp(i,j,:) = obc_temp;
        salt(i,j,:) = obc_salt;
    end
end
netcdf.putVar(nc,obc_temp_varid,temp);
netcdf.putVar(nc,obc_salinity_varid,salt);

% close file
netcdf.close(nc);

if(ftbverbose); fprintf(['end   : ' subname '\n']);end;
