function write_FVCOM_elevtide(Mobj,MJD,ElevationFile,MyTitle)
% Write an FVCOM surface elevation time series forcing file 
%
% write_FVCOM_elevtide(Mobj,MJD,ElevationFile,MyTitle)
%
% DESCRIPTION:
%    Write an FVCOM NetCDF surface elevation forcing file
%
% INPUT:
%   Mobj         = Matlab mesh object.
%   MJD          = list of modified Modified Julian Dates of size [times]
%                   (defined as unlimited in the NetCDF file).
%   ElevationFile    = name of NetCDF file.
%   MyTitle      = casename title, written as global attribute of NetCDF file.
%
% OUTPUT:
%    ElevationFile, A NetCDF FVCOM surface elevations tide forcing file
%
% EXAMPLE USAGE
%    write_FVCOM_elevtide(Mobj,MJD,ElevationFile,MyTitle)
%
% Author(s):  
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Karen Thurston (National Oceanography Centre Liverpool)
% 
% Revision history
%    2012-08-08 (PWC) First version.
%    2012-11-14 (PWC) Updated to expect Modified Julian Day rather than doing 
%    the conversion in here. Also put the pieces in set_elevtide in here to
%    simplify the process of writing out an elevation input file.
%    2012-12-04 (KJT) Updated to use surface elevation and open boundary 
%    nodes from Mobj.
%   
%==============================================================================

global ftbverbose 
report = false;
if(ftbverbose); report = true; end
subname = 'write_FVCOM_elevtide';
if(report); fprintf('\n'); end
if(report); fprintf(['begin : ' subname '\n']); end

% Get a list of the open boundary nodes. Transpose Mobj.obc_nodes so the
% order of the boundary nodes is preserved.
tmpObcNodes = Mobj.obc_nodes';
% Flip it back so it's the same shape as it would have been using the old
% code.
ObcNodes = tmpObcNodes(tmpObcNodes~=0)';

%------------------------------------------------------------------------------
% Sanity check on input and dimensions
%------------------------------------------------------------------------------
nTimes = numel(MJD);
if(report); fprintf('Number of time steps %d\n',nTimes); end

nObcs = numel(ObcNodes);
if(report); fprintf('Number of Open Boundary Nodes %d\n',nObcs); end

[chk1, chk2] = size(Mobj.surfaceElevation);
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
date_str_len_dimid=netcdf.defDim(nc,'DateStrLen',26);

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
netcdf.putVar(nc,time_varid,0,nTimes,MJD);
netcdf.putVar(nc,itime_varid,floor(MJD));
netcdf.putVar(nc,itime2_varid,0,nTimes,mod(MJD,1)*24*3600*1000);
nStringOut = char();
for i=1:nTimes
    [nYr, nMon, nDay, nHour, nMin, nSec] = mjulian2greg(MJD(i));
    if strcmp(sprintf('%02i', nSec), '60')
        % Fix some weirdness with mjulian2greg. I think this is caused by
        % rounding errors. My testing suggests this is not a problem around
        % midnight, so the number of days (and thus possibly months and
        % years) is unaffected.
        if mod(nMin + 1, 60) == 0
            % Up the hour by one too
            nHour = mod(nHour + 1, 24);
        end
        nMin = mod(nMin + 1, 60);
        nSec = 0;
    end
    nDate = [nYr, nMon, nDay, nHour, nMin, nSec];
    nStringOut = [nStringOut, sprintf('%04i/%02i/%02i %02i:%02i:%02i       ',nDate)];
end
netcdf.putVar(nc,Times_varid,nStringOut);
netcdf.putVar(nc,elevation_varid,Mobj.surfaceElevation);

% close file
netcdf.close(nc);

if(report); fprintf(['end   : ' subname '\n']); end;

