function write_FVCOM_river_ERSEM(RiverFile,RiverName,time,flux,temp,salt,n1p,n3n,n4n,n5s,dic,o2,bioalk,alkalinity,RiverInfo1,RiverInfo2)
% Write FVCOM 3.x netCDF river file
%
% function write_FVCOM_river_ERSEM(RiverFile,RiverName,time,flux,temp,salt,n1p,n3n,n4n,n5s,dic,alkalinity,RiverInfo1,RiverInfo2)
%
% DESCRIPTION:
%    Write river flux, temperature, salinity and ERSEM nutrients to an
%    FVCOM river file.
%
% INPUT
%    RiverFile  : FVCOM 3.x netCDF river forcing file
%    RiverName  : Name of the actual River
%    time       : Timestamp array in modified Julian day
%    flux       : Total river flux in m^3/s (dimensions [time, nRivernodes])
%    temp       : Temperature in C (dimensions [time, nRivernodes])
%    salt       : Salinity in PSU (dimensions [time, nRivernodes])
%    n1p        : Phosphate (mmol P/m^3) (dimensions [time, nRivernodes])
%    n3n        : Nitrate (mmol N/m^3) (dimensions [time, nRivernodes])
%    n4n        : Ammonia (mmol N/m^3) (dimensions [time, nRivernodes])
%    n5s        : Silicate (mmol Si/m^3) (dimensions [time, nRivernodes])
%    dic        : Dissolved inorganic carbon (mmol C/m^3) (dimensions [time, nRivernodes])
%    O2        : Dissolved Oxygen (mmol O_2/m^3) (dimensions [time, nRivernodes])
%    bioalk        : Bio_alkalinity compartment carbon (mmol C/m^3) (dimensions [time, nRivernodes])
%    Total_alk : Total Alkalinity (umol C/m^3) (dimensions [time, nRivernodes])
%    RiverInfo1 : Global attribute title of file
%    RiverInfo2 : Global attribute info of file
%
% OUTPUT:
%    FVCOM netCDF river file with flux, temperature, salinity and ERSEM
%    nutrients.
%
% EXAMPLE USAGE
%    write_FVCOM_river('tst_riv.nc', {'Penobscot'}, time, flux, ...
%         temp, salt, n1p, n3n, n4n, n5s, dic, alkalinity, ...
%         'Penobscot Flux', 'source: USGS')
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Ricardo Torres (Plymouth Marine Laboratory)
%
% Revision history
%   2016-03-14 New version to export nutrients alongside the physical
%   parameters for FVCOM-ERSEM. Based on write_FVCOM_river.
%
%==========================================================================

[~, subname] = fileparts(mfilename('fullpath'));

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

[nTimes, nRivnodes] = size(flux);
if ftbverbose
    fprintf('# of river nodes: %d\n', nRivnodes);
    fprintf('# of time frames: %d\n', nTimes);
end

[year1, month1, day1, ~, ~, ~] = mjulian2greg(time(1));
[year2, month2, day2, ~, ~, ~] = mjulian2greg(time(end));
if ftbverbose
    fprintf('time series begins at:\t%04d %02d %02d\n', year1, month1, day1);
    fprintf('time series ends at:\t%04d %02d %02d\n', year2, month2, day2);
end
clear year? month? day?

%--------------------------------------------------------------------------
% dump to netcdf file
%--------------------------------------------------------------------------

% river node forcing
nc = netcdf.create(RiverFile, 'clobber');

% global variables
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'type', 'FVCOM RIVER FORCING FILE')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'title', RiverInfo1)
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'info', RiverInfo2)
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'history', sprintf('File created using %s from the MATLAB fvcom-toolbox', subname))

% dimensions
namelen_dimid = netcdf.defDim(nc, 'namelen', 80);
rivers_dimid = netcdf.defDim(nc, 'rivers', nRivnodes);
time_dimid = netcdf.defDim(nc, 'time', netcdf.getConstant('NC_UNLIMITED'));
date_str_len_dimid = netcdf.defDim(nc, 'DateStrLen', 26);

% variables
river_names_varid = netcdf.defVar(nc, 'river_names', 'NC_CHAR', [namelen_dimid, rivers_dimid]);

time_varid = netcdf.defVar(nc, 'time', 'NC_FLOAT', time_dimid);
netcdf.putAtt(nc, time_varid, 'long_name', 'time');
netcdf.putAtt(nc, time_varid, 'units', 'days since 1858-11-17 00:00:00');
netcdf.putAtt(nc, time_varid, 'format', 'modified julian day (MJD)');
netcdf.putAtt(nc, time_varid, 'time_zone', 'UTC');

itime_varid = netcdf.defVar(nc, 'Itime', 'NC_INT', time_dimid);
netcdf.putAtt(nc, itime_varid, 'units', 'days since 1858-11-17 00:00:00');
netcdf.putAtt(nc, itime_varid, 'format', 'modified julian day (MJD)');
netcdf.putAtt(nc, itime_varid, 'time_zone', 'UTC');

itime2_varid = netcdf.defVar(nc, 'Itime2', 'NC_INT', time_dimid);
netcdf.putAtt(nc, itime2_varid, 'units', 'msec since 00:00:00');
netcdf.putAtt(nc, itime2_varid, 'time_zone', 'UTC');

times_varid = netcdf.defVar(nc,'Times','NC_CHAR',[date_str_len_dimid, time_dimid]);
netcdf.putAtt(nc, times_varid, 'time_zone','UTC');

river_flux_varid = netcdf.defVar(nc, 'river_flux', 'NC_FLOAT', [rivers_dimid, time_dimid]);
netcdf.putAtt(nc, river_flux_varid, 'long_name', 'river runoff volume flux');
netcdf.putAtt(nc, river_flux_varid, 'units', 'm^3s^-1');

river_temp_varid = netcdf.defVar(nc, 'river_temp', 'NC_FLOAT', [rivers_dimid, time_dimid]);
netcdf.putAtt(nc, river_temp_varid, 'long_name', 'river runoff temperature');
netcdf.putAtt(nc, river_temp_varid, 'units', 'Celsius');

river_salt_varid = netcdf.defVar(nc, 'river_salt', 'NC_FLOAT', [rivers_dimid, time_dimid]);
netcdf.putAtt(nc, river_salt_varid, 'long_name', 'river runoff salinity');
netcdf.putAtt(nc, river_salt_varid, 'units', 'PSU');

river_n1p_varid = netcdf.defVar(nc, 'N1_p', 'NC_FLOAT', [rivers_dimid, time_dimid]);
netcdf.putAtt(nc, river_n1p_varid, 'long_name', 'phosphate phosphorus');
netcdf.putAtt(nc, river_n1p_varid, 'units', 'mmol P/m^3');

river_n3n_varid = netcdf.defVar(nc, 'N3_n', 'NC_FLOAT', [rivers_dimid, time_dimid]);
netcdf.putAtt(nc, river_n3n_varid, 'long_name', 'nitrate nitrogen');
netcdf.putAtt(nc, river_n3n_varid, 'units', 'mmol N/m^3');

river_n4n_varid = netcdf.defVar(nc, 'N4_n', 'NC_FLOAT', [rivers_dimid, time_dimid]);
netcdf.putAtt(nc, river_n4n_varid, 'long_name', 'ammonium nitrogen');
netcdf.putAtt(nc, river_n4n_varid, 'units', 'mmol N/m^3');

river_n5s_varid = netcdf.defVar(nc, 'N5_s', 'NC_FLOAT', [rivers_dimid, time_dimid]);
netcdf.putAtt(nc, river_n5s_varid, 'long_name', 'silicate silicate');
netcdf.putAtt(nc, river_n5s_varid, 'units', 'mmol Si/m^3');

river_O2_varid = netcdf.defVar(nc, 'O2_o', 'NC_FLOAT', [rivers_dimid, time_dimid]);
netcdf.putAtt(nc, river_O2_varid, 'long_name', 'dissolved Oxygen');
netcdf.putAtt(nc, river_O2_varid, 'units', 'mmol O_2/m^3');

river_TA_varid = netcdf.defVar(nc, 'O3_TA', 'NC_FLOAT', [rivers_dimid, time_dimid]);
netcdf.putAtt(nc, river_TA_varid, 'long_name', 'carbonate total alkalinity');
netcdf.putAtt(nc, river_TA_varid, 'units', 'mmol C/m^3');


river_dic_varid = netcdf.defVar(nc, 'O3_c', 'NC_FLOAT', [rivers_dimid, time_dimid]);
netcdf.putAtt(nc, river_dic_varid, 'long_name', 'carbonate total dissolved inorganic carbon');
netcdf.putAtt(nc, river_dic_varid, 'units', 'mmol C/m^3');

river_bioalk_varid = netcdf.defVar(nc, 'O3_bioalk', 'NC_FLOAT', [rivers_dimid, time_dimid]);
netcdf.putAtt(nc, river_bioalk_varid, 'long_name', 'carbonate bioalkalinity');
netcdf.putAtt(nc, river_bioalk_varid, 'units', 'umol/kg');

% end definitions
netcdf.endDef(nc);

% river names (must be 80 character strings)
rString = char();
for i = 1:nRivnodes
    % Left-aligned 80 character string.
    rString = [rString, sprintf('%-80s', RiverName{i})];
end
netcdf.putVar(nc, river_names_varid, rString);

% dump dynamic data
netcdf.putVar(nc, time_varid, 0, nTimes, time);
netcdf.putVar(nc, itime_varid, 0, nTimes, floor(time));
netcdf.putVar(nc, itime2_varid, 0, nTimes, mod(time, 1)*24*3600*1000);
netcdf.putVar(nc, river_flux_varid, flux');
netcdf.putVar(nc, river_temp_varid, temp');
netcdf.putVar(nc, river_salt_varid, salt');
netcdf.putVar(nc, river_n1p_varid, n1p');
netcdf.putVar(nc, river_n3n_varid, n3n');
netcdf.putVar(nc, river_n4n_varid, n4n');
netcdf.putVar(nc, river_n5s_varid, n5s');
netcdf.putVar(nc, river_dic_varid, dic');
netcdf.putVar(nc, river_O2_varid, o2');
netcdf.putVar(nc, river_TA_varid, alkalinity');
netcdf.putVar(nc, river_bioalk_varid,bioalk');
% build the time string and output to netCDF.
nStringOut = char();
[nYr, nMon, nDay, nHour, nMin, nSec] = mjulian2greg(time);
for tt = 1:nTimes
    nDate = [nYr(tt), nMon(tt), nDay(tt), nHour(tt), nMin(tt), nSec(tt)];
    nStringOut = [nStringOut, sprintf('%04i/%02i/%02i %02i:%02i:%02i       ', nDate)];
end
netcdf.putVar(nc, times_varid, nStringOut);

netcdf.close(nc);

if ftbverbose
    fprintf('end   : %s\n', subname)
end

