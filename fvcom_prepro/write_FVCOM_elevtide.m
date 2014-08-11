function write_FVCOM_elevtide(Mobj,MJD,ElevationFile,MyTitle,varargin)
% Write an FVCOM surface elevation time series forcing file 
%
% write_FVCOM_elevtide(Mobj, MJD, ElevationFile, MyTitle)
%
% DESCRIPTION:
%    Write an FVCOM NetCDF surface elevation forcing file
%
% INPUT:
%   Mobj = Matlab mesh object with fields:
%       obc_nodes - array of boundary node IDs.
%       surfaceElevation - array of surface elevation values (shaped [space,
%       time]).
%   MJD = list of modified Modified Julian Dates of size [times] (defined
%       as unlimited in the netCDF file).
%   ElevationFile = name of netCDF file.
%   MyTitle = casename title, written as global attribute of netCDF file.
%   Optional keyword-argument pairs. These control the time variables. This
%   script defaults to writing 'Times' only.
%   FVCOM needs only one of:
%       1. Times: character string of times
%       2. Itime and Itime2: integer days and milliseconds since midnight
%       3. time: float days.
%   FVCOM checks for these in the order above and this script defaults to
%   writing Times only. Adjust the keyword-argument pairs to your liking:
%
%   'strtime' = set to true to output the 'Times' variable
%   'inttime' = set to true to output the 'Itime' and 'Itime2' variables
%   'floattime' = set to true to output the 'time' variable
%
% OUTPUT:
%    ElevationFile, A NetCDF FVCOM surface elevations tide forcing file
%
% EXAMPLE USAGE
%   With default settings:
%       write_FVCOM_elevtide(Mobj, MJD, '/tmp/elevtide.nc, 'Shelf tides')
%   Enable the 'time' variable in the netCDF.
%       write_FVCOM_elevtide(Mobj, MJD, '/tmp/elevtide.nc, ...
%           'Shelf tides', 'floattime', true)
%
% Author(s):  
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Karen Amoudry (National Oceanography Centre Liverpool)
% 
% Revision history
%   2012-08-08 (PWC) First version.
%   2012-11-14 (PWC) Updated to expect Modified Julian Day rather than
%   doing the conversion in here. Also put the pieces in set_elevtide in
%   here to simplify the process of writing out an elevation input file.
%   2012-12-04 (KJA) Updated to use surface elevation and open boundary
%   nodes from Mobj.
%   2013-08-16 (KJA) Updated output of Itime2 to avoid rounding errors when
%   converting from double to single format.
%   2013-09-03 - Removed PWC's fix for timestrings. Issue was due to
%   rounding errors caused by mjulian2greg.m, which have now been fixed.
%   2014-01-27 - (PWC) Simplify the ftbverbose/report stuff.
%   2014-08-11 - (PWC) Add new flags to control which time variables to
%   use. FVCOM reads the 'Times' variable first if present, then falls back
%   to 'Itime' and 'Itime2' and finally 'time'. Also reinstate the original
%   version of the calculation of Itime2 as the rounding effect was
%   smoothing out the data too much, affecting its precision.
%   
%==========================================================================

global ftbverbose

subname = 'write_FVCOM_elevtide';
if ftbverbose; fprintf('\nbegin : %s \n', subname); end

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

% Get a list of the open boundary nodes. Transpose Mobj.obc_nodes so the
% order of the boundary nodes is preserved.
tmpObcNodes = Mobj.obc_nodes';
% Flip it back so it's the same shape as it would have been using the old
% code.
ObcNodes = tmpObcNodes(tmpObcNodes~=0)';

%--------------------------------------------------------------------------
% Sanity check on input and dimensions
%--------------------------------------------------------------------------
nTimes = numel(MJD);
if ftbverbose; fprintf('Number of time steps %d\n',nTimes); end

nObcs = numel(ObcNodes);
if ftbverbose; fprintf('Number of Open Boundary Nodes %d\n',nObcs); end

[chk1, chk2] = size(Mobj.surfaceElevation);
if nObcs ~= chk1 || nTimes ~= chk2
    fprintf('Surface elevation dimensions do not match time series and number of boundary nodes.\n')
    fprintf('Surface elevation nodes and time sizes: (%d, %d)\n', chk1, chk2)
    fprintf('Boundary nodes size: %d\n', nObcs)
    fprintf('Times size: %d\n', nTimes)
	error('Input data sizes do not match. Check and try again.');
end

%%
%--------------------------------------------------------------------------
% Dump the file
%--------------------------------------------------------------------------

nc=netcdf.create(ElevationFile,'clobber');

% define global attributes
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type','FVCOM TIME SERIES ELEVATION FORCING FILE')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title',MyTitle)
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history','File created using write_FVCOM_elevtide from the MATLAB fvcom-toolbox')

% define dimensions
nobc_dimid=netcdf.defDim(nc,'nobc',nObcs);
time_dimid=netcdf.defDim(nc,'time',netcdf.getConstant('NC_UNLIMITED'));
date_str_len_dimid=netcdf.defDim(nc,'DateStrLen',26);

% define variables and attributes
nobc_varid=netcdf.defVar(nc,'obc_nodes','NC_INT',nobc_dimid);
netcdf.putAtt(nc,nobc_varid,'long_name','Open Boundary Node Number');
netcdf.putAtt(nc,nobc_varid,'grid','obc_grid');

iint_varid=netcdf.defVar(nc,'iint','NC_INT',time_dimid);
netcdf.putAtt(nc,iint_varid,'long_name','internal mode iteration number');

if floattime
    time_varid=netcdf.defVar(nc,'time','NC_FLOAT',time_dimid);
    netcdf.putAtt(nc,time_varid,'long_name','time');
    netcdf.putAtt(nc,time_varid,'units','days since 1858-11-17 00:00:00');
    netcdf.putAtt(nc,time_varid,'format','modified julian day (MJD)');
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

if strtime
    Times_varid=netcdf.defVar(nc,'Times','NC_CHAR',[date_str_len_dimid, time_dimid]);
    netcdf.putAtt(nc,Times_varid,'time_zone','UTC');
end

elevation_varid=netcdf.defVar(nc,'elevation','NC_FLOAT',[nobc_dimid, time_dimid]);
netcdf.putAtt(nc,elevation_varid,'long_name','Open Boundary Elevation');
netcdf.putAtt(nc,elevation_varid,'units','meters');

% end definitions
netcdf.endDef(nc);

% write data
netcdf.putVar(nc,nobc_varid,ObcNodes);
netcdf.putVar(nc,iint_varid,0,nTimes,1:nTimes);
if floattime
    netcdf.putVar(nc,time_varid,0,nTimes,MJD);
end
if inttime
    netcdf.putVar(nc,itime_varid,floor(MJD));
    netcdf.putVar(nc,itime2_varid,0,nTimes,round(mod(MJD,1) * 24 * 60 * 60 * 1000));
end
if strtime
    nStringOut = char();
    [nYr, nMon, nDay, nHour, nMin, nSec] = mjulian2greg(MJD);
    for i=1:nTimes
        nDate = [nYr(i), nMon(i), nDay(i), nHour(i), nMin(i), nSec(i)];
        nStringOut = [nStringOut, sprintf('%04i/%02i/%02i %02i:%02i:%09.6f', nDate)];
    end
    netcdf.putVar(nc,Times_varid,nStringOut);
end
netcdf.putVar(nc,elevation_varid,Mobj.surfaceElevation);

% close file
netcdf.close(nc);

if ftbverbose; fprintf('end   : %s \n', subname); end

