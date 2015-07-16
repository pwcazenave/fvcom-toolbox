function write_FVCOM_tsobc(basename,time,nSiglay,in_temp,in_salt,Mobj,varargin)
% Export temperature and salinity forcing at the open boundary.
%
% function write_FVCOM_tsobc(basename,time,nSiglay,in_temp,in_salt)
%
% DESCRIPTION:
%    Setup an FVCOM hydrographic open boundary forcing file. Supply either
%    uniform values for temperature and salinity or 3D arrays (node,
%    sigma_layers, time).
%
% INPUT
%    basename - Model case name (to find the bathymetry and open boundary
%    .dat files).
%    time - Time (Modified Julian Days)
%    nSiglay - Number of sigma layers
%    in_temp - Boundary temperature (Celsius)
%    in_salt - Boundary salinity (psu)
%    Mobj - Mesh Object
%    Optional keyword-argument pairs. These control the time variables.
%    This script defaults to writing 'Times' only.
%    FVCOM needs only one of:
%        1. Times: character string of times
%        2. Itime and Itime2: integer days and milliseconds since midnight
%        3. time: float days.
%    FVCOM checks for these in the order above and this script defaults to
%    writing Times only. Adjust the keyword-argument pairs to your liking:
%
%    'strtime' = set to true to output the 'Times' variable
%    'inttime' = set to true to output the 'Itime' and 'Itime2' variables
%    'floattime' = set to true to output the 'time' variable
%
% OUTPUT:
%    FVCOM hydrographic open boundary netCDF file
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Karen Amoudry (National Oceanography Centre, Liverpool)
%
% PWC Revision history
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
%
% KJA Revision history
%    Undated - Add better check for the size of the input arrays (works
%    with scalars).
%    2013-08-16 - Updated output of Itime2 to avoid rounding errors
%    when converting from double to single format.
%
%==============================================================================

if nargin == 5
    warning(['Assuming uniform terrain-following sigma coordinates. ', ...
        'To specify different sigma coordintes, supply a MATLAB mesh ', ...
        'structure with fields siglay and siglev.'])
end

subname = 'write_FVCOM_tsobc';
global ftbverbose;
if ftbverbose
  fprintf('\nbegin : %s\n', subname)
end

% Default to string times as FVCOM looks for these first.
strtime = true;
inttime = false;
floattime = false;
for vv = 1:2:length(varargin)
    switch varargin{vv}
        case 'strtime'
            strtime = true;
        case 'inttime'
            inttime = true;
        case 'floattime'
            floattime = true;
    end
end

fvcom_bathy = [basename, '_dep.dat'];
fvcom_obc   = [basename, '_obc.dat'];
tsOBCFile = [basename, '_tsobc.nc'];

%------------------------------------------------------------------------------
% read in the FVCOM open boundary node data (need node numbers and dimension)
%------------------------------------------------------------------------------
fid = fopen(fvcom_obc,'r');
if(fid  < 0)
  error(['file: ' fvcom_obc ' does not exist']);
end
C = textscan(fid, '%s %s %s %s %d', 1);
nObc = C{5};
obc_nodes = zeros(nObc,1);
if(ftbverbose); fprintf('reading obc file\n'); end;
if(ftbverbose); fprintf('# nodes %d\n',nObc); end;
for i=1:nObc
  C = textscan(fid, '%d %d %d', 1);
  obc_nodes(i) = C{2};
end

if(ftbverbose); fprintf('obc reading complete\n');end;

%------------------------------------------------------------------------------
% read in the FVCOM bathymetry data (need bathymetry on open boundary nodes)
%------------------------------------------------------------------------------
fid = fopen(fvcom_bathy,'r');
if(fid  < 0)
  error(['file: ' fvcom_bathy ' does not exist']);
end
C = textscan(fid, '%s %s %s %d', 1);
Nverts = C{4};
h = zeros(Nverts,1);
if(ftbverbose); fprintf('reading bathymetry file\n');end;
if(ftbverbose); fprintf('# nodes %d\n',Nverts);end;
for i=1:Nverts
  C = textscan(fid, '%f %f %f', 1);
  h(i) = C{3};
end
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
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title','Open boundary temperature and salinity nudging')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type','FVCOM TIME SERIES OBC TS FILE')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history','File generated using write_FVCOM_tsobc.m from the MATLAB fvcom-toolbox')


% define dimensions
nobc_dimid=netcdf.defDim(nc,'nobc',nObc);
datestrlen_dimid=netcdf.defDim(nc,'DateStrLen',26);
time_dimid=netcdf.defDim(nc,'time',netcdf.getConstant('NC_UNLIMITED'));
siglay_dimid=netcdf.defDim(nc,'siglay',nSiglay);
siglev_dimid=netcdf.defDim(nc,'siglev',nSiglev);

% variables
if strtime
    Times_varid=netcdf.defVar(nc,'Times','NC_CHAR',[datestrlen_dimid, time_dimid]);
    netcdf.putAtt(nc,Times_varid,'time_zone','UTC');
end

if floattime
    time_varid=netcdf.defVar(nc,'time','NC_FLOAT',time_dimid);
    netcdf.putAtt(nc,time_varid,'long_name','time');
    netcdf.putAtt(nc,time_varid,'units','days since 1858-11-17 00:00:00');
    netcdf.putAtt(nc,time_varid,'time_zone','UTC');
end

if inttime
    itime_varid=netcdf.defVar(nc,'Itime','NC_INT',time_dimid);
    netcdf.putAtt(nc,itime_varid,'units','days since 1858-11-17 00:00:00');
    netcdf.putAtt(nc,itime_varid,'format','modified julian day (MJD)');
    netcdf.putAtt(nc,itime_varid,'time_zone','UTC');

    itime2_varid=netcdf.defVar(nc,'Itime2','NC_INT',time_dimid);
    netcdf.putAtt(nc,itime2_varid,'units','msec since 00:00:00');
    netcdf.putAtt(nc,itime2_varid,'time_zone','UTC');
end

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
if strtime
    nStringOut = char();
    [nYr, nMon, nDay, nHour, nMin, nSec] = mjulian2greg(time);
    for i=1:nTimes
        nDate = [nYr(i), nMon(i), nDay(i), nHour(i), nMin(i), nSec(i)];
        nStringOut = [nStringOut, sprintf('%04i/%02i/%02i %02i:%02i:%09.6f', nDate)];
    end
    netcdf.putVar(nc,Times_varid,[0, 0], [26, nTimes],nStringOut);
end
if floattime
    netcdf.putVar(nc,time_varid,0,numel(time),time);
end
if inttime
    netcdf.putVar(nc,itime_varid,floor(time));
    %netcdf.putVar(nc,itime2_varid,0,numel(time),mod(time,1)*24*3600*1000); % PWC original
    % KJA edit: avoids rounding errors when converting from double to single
    % Rounds to nearest multiple of the number of msecs in an hour
    netcdf.putVar(nc,itime2_varid,0,numel(time),round((mod(time,1)*24*3600*1000)/(3600*1000))*(3600*1000));
end
netcdf.putVar(nc,obc_temp_varid,temp);
netcdf.putVar(nc,obc_salinity_varid,salt);

% close file
netcdf.close(nc);

if ftbverbose; fprintf('end   : %s\n', subname); end
