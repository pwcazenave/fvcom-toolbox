function write_FVCOM_tsobc(basename,time,nSiglay,in_temp,in_salt,Mobj)
% example file for dumping a file to force temperature and salinity at the open b.
%
% function write_FVCOM_tsobc(basename,time,nSiglay,in_temp,in_salt)
%
% DESCRIPTION:
%    Setup an FVCOM hydrographic open boundary forcing file. Supply either
%    uniform values for temperature and salinity or 3D arrays (node,
%    sigma_layers, time).
%
% INPUT
%    Model case name (to find the bathymetry and open boundary .dat files).
%    Time
%    Number of sigma layers
%    Boundary temperature (Celcius)
%    Boundary salinity (psu)
%    Mobj (optional)
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
%    2012-10-08 Updated help to reflect the fact nSiglev is calculated as
%    nSiglay+1.
%    2012-11-09 Add new arguments to use user defined temperature and
%    salinity.
%    2013-01-09 Add support for 3D input temperature and salinity (such as
%    might be generated with get_POLCOMS_tsobc.m.
%    KJT: Add better check for the size of the input arrays (works with 
%    scalars).
%
%==============================================================================

if nargin == 5
    warning(['Assuming uniform terrain-following sigma coordinates. ',...
        'To specify different sigma coordintes, supply a MATLAB mesh ',...
        'structure with fields siglay and siglev.'])
end

subname = 'write_FVCOM_tsobc';
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
% Generate the requisite data
%--------------------------------------------------------------

% extract bathymetry at open boundary nodes
obc_h = h(obc_nodes);

% time
% time = 0:1:31.;
nTimes = numel(time);

nSiglev = nSiglay + 1;

% Create or process the temperature and salinity arrays.
if max(size(in_temp)) == 1
    inc = 1/real(nSiglay);
    siglev = 0:-inc:-1;
    siglay = nan(1, nSiglay);
    for i=1:nSiglay
        siglay(i) = mean(siglev(i:i+1));
    end
    % initialize temperature/salinity arrays
    temp = zeros(nObc,nSiglay,nTimes);
    salt = zeros(nObc,nSiglay,nTimes);

    % set a constant temperature and salinity
    obc_temp = repmat(in_temp, 1, nTimes);
    obc_salt = repmat(in_salt, 1, nTimes);

    % set variable temperature and salinity
    % for i=1:nTimes
    % 	obc_temp(i) = 18. + 2.*real(i-1)/nTimes;
    % 	obc_salt(i) = 30. - 5.*real(i-1)/nTimes;
    % end

    % Create 3D array from three 1D arrays
    % temp = repmat(obc_temp, [nObc, nSiglay, 1]);
    % salt = repmat(obc_salt, [nObc, nSiglay, 1]);
    for i=1:nObc
        for j=1:nSiglay
            temp(i,j,:) = obc_temp;
            salt(i,j,:) = obc_salt;
        end
    end
else
    % We have a 3D array already so we just need a couple of stats.
    temp = in_temp;
    salt = in_salt;

    if nargin == 6 && isfield(Mobj, 'siglay') && isfield(Mobj, 'siglev')
        siglev = Mobj.siglev;
        siglay = Mobj.siglay;
    else
        warning('Assuming uniform terrain-following sigma coordinates')
        inc = 1/real(nSiglay);
        siglev = 0:-inc:-1;
        siglay = nan(1, nSiglay);
    end

    if nSiglev ~= size(in_temp, 2) + 1 || length(siglev) ~= size(in_temp, 2) + 1 || length(siglev) ~= size(in_salt, 2) + 1
        error('Specified number sigma levels does not match supplied data')
    end
    if nSiglay ~= size(in_temp, 2) || length(siglay) ~= size(in_temp, 2) || length(siglay) ~= size(in_salt, 2)
        error('Specified number of sigma layers does not match supplied data')
    end
end

%--------------------------------------------------------------
% set NetCDF variables and dump to file
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
netcdf.putAtt(nc,time_varid,'units','days since 1858-11-17 00:00:00');
netcdf.putAtt(nc,time_varid,'time_zone','UTC');

itime_varid=netcdf.defVar(nc,'Itime','NC_INT',time_dimid);
netcdf.putAtt(nc,itime_varid,'units','days since 1858-11-17 00:00:00');
netcdf.putAtt(nc,itime_varid,'format','modified julian day (MJD)');
netcdf.putAtt(nc,itime_varid,'time_zone','UTC');

itime2_varid=netcdf.defVar(nc,'Itime2','NC_INT',time_dimid);
netcdf.putAtt(nc,itime2_varid,'units','msec since 00:00:00');
netcdf.putAtt(nc,itime2_varid,'time_zone','UTC');

nobc_varid=netcdf.defVar(nc,'obc_nodes','NC_INT',nobc_dimid);
netcdf.putAtt(nc,nobc_varid,'long_name','Open Boundary Node Number');
netcdf.putAtt(nc,nobc_varid,'grid','obc_grid');
netcdf.putAtt(nc,nobc_varid,'type','data');

obc_h_varid=netcdf.defVar(nc,'obc_h','NC_FLOAT',nobc_dimid);
netcdf.putAtt(nc,obc_h_varid,'long_name','Open Boundary Depth');
netcdf.putAtt(nc,obc_h_varid,'units','m');
netcdf.putAtt(nc,obc_h_varid,'grid','obc_grid');
netcdf.putAtt(nc,obc_h_varid,'type','data');

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

netcdf.putVar(nc,obc_temp_varid,temp);
netcdf.putVar(nc,obc_salinity_varid,salt);

% close file
netcdf.close(nc);

if(ftbverbose); fprintf(['end   : ' subname '\n']);end;
