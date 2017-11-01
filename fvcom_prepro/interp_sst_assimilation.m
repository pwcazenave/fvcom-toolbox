function Mobj = interp_sst_assimilation(Mobj, conf, output_file)
% Interpolate SST values to a given FVCOM mesh object.
%
% Mobj = interp_sst_assimilation(Mobj, year, sst_dir, file_pattern)
%
% DESCRIPTION:
%   Interpolate SST data from remote sensing data onto the supplied model
%   grid.
%
% INPUT:
%   Mobj - MATLAB mesh object containing fields:
%       lon, lat - node coordinates (spherical)
%       lonc, latc - element coordinates (spherical)
%       tri - element triangulation table
%   conf - struct with fields:
%       sst_dir - directory containing the SST data
%       sst_pattern - file name pattern for the SST data
%       year - year for which to generate SST data
%   output_file - path to which to output the netCDF
%
% OUTPUT:
%   FVCOM data assimilation SST netCDF file.
%   Mobj - input MATLAB mesh object with added 'assim.sst' field, with
%   fields:
%       data - the SST data
%       time - the SST time series
%
% EXAMPLE USAGE:
%   Mobj = read_fvcom_mesh('casename.grd');
%   conf.sst_dir = '/home/user/GHRSST/';
%   conf.sst_pattern = '-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc';
%   interp_sst_assimilation(Mobj, conf, 'casename_sstgrd.nc');
%
% Author(s)
%   Ricardo Torres (Plymouth Marine Laboratory)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% History:
%   2017-08-02 Turned the script into a function.

[~, subname] = fileparts(mfilename('fullpath'));
global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

year = conf.year;
% SST files.
sst_dir = conf.sst_dir;
sst_pattern = conf.sst_pattern;

time_span = datenum(year-1,12,31):datenum(year+1,01,01);
file_list = cell(length(time_span), 1);
% Add last day from the previous year and the first day of the following
% year.
for aa = 1:length(time_span)
    file_list{aa} = fullfile(sst_dir, ...
        datestr(time_span(aa),'yyyy'), ...
        sprintf(sst_pattern, ...
        str2num(datestr(time_span(aa), 'yyyy')), ...
        str2num(datestr(time_span(aa), 'mm')), ...
        str2num(datestr(time_span(aa), 'dd'))));
    if ~exist(file_list{aa}, 'file')
        error('We are missing a file (%s) from this year (%04d)', ...
            file_list{aa}, datestr(time_span(aa),'yyyy'))
    end
end

% Read SST data files and interpolate each to the FVCOM mesh
lon = ncread(file_list{1},'lon');
lat = ncread(file_list{1},'lat');
mask = ncread(file_list{1},'mask');
[lonm,latm]=meshgrid(lon,lat);
lonm=lonm';
latm=latm';
lonm = lonm(mask==1);
latm = latm(mask==1);
time = zeros(length(file_list),1);
sst = zeros(Mobj.nVerts,length(file_list),1,'single');
fvcomlon = Mobj.lon;
fvcomlat = Mobj.lat;

if license('test', 'Distrib_Computing_Toolbox')
    if isempty(gcp('nocreate'))
        % Force pool to be local in case we have remote pools available.
        parpool('local');
    end
end

if ftbverbose
    fprintf('Progress:\n');
    fprintf([repmat('.', 1, length(file_list)), '\n']);
end
parfor ff = 1:length(file_list)
    if ftbverbose
        fprintf('|');
    end
    sst_eo = ncread(file_list{ff}, 'analysed_sst') - 273.15;
    mask = ncread(file_list{ff}, 'mask');
    lon = ncread(file_list{ff}, 'lon');
    lat = ncread(file_list{ff}, 'lat');
    [lonm, latm] = meshgrid(lon, lat);
    lonm = lonm';
    latm = latm';
    lonm = lonm(mask == 1);
    latm = latm(mask == 1);

    time_eo = ncread(file_list{ff}, 'time');
    time_eo_units = ncreadatt(file_list{ff}, 'time', 'units');
    t0str = textscan(time_eo_units, 'seconds since %s%s');

    t0 = datenum([strtrim(t0str{1}{1}), strtrim(t0str{2}{1})], 'yyyy-mm-ddHH:MM:SS');
    time_out = (t0 + double(time_eo/(60*60*24)));
    sst_eo = sst_eo(mask == 1);

    % Build interpolant
    ft = scatteredInterpolant(double(lonm), double(latm), sst_eo, 'nearest', 'linear');

    sst(:,ff) = ft(fvcomlon, fvcomlat);
    time(ff) = time_out + 0.5; % fvcom expects these to be at mid-day
end

ntimes = length(time);

% Do the times.
[sYr, sMon, sDay, sHr, sMin, sSec] = datevec(time);
MJDtime = greg2mjulian(sYr, sMon, sDay, sHr, sMin, sSec);

% Create netCDF file
nc = netcdf.create(output_file, 'clobber');

netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'year', num2str(year))
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title','FVCOM SST 1km merged product File')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'institution','Plymouth Marine Laboratory')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'source','FVCOM grid (unstructured) surface forcing')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history', sprintf('File created with %s from the MATLAB fvcom-toolbox', subname))
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'references','http://fvcom.smast.umassd.edu, http://codfish.smast.umassd.edu')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.0')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'CoordinateSystem',Mobj.nativeCoords)
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'CoordinateProjection','init=WGS84') % WGS84?

% Dimensions
nele_dimid=netcdf.defDim(nc,'nele',Mobj.nElems);
node_dimid=netcdf.defDim(nc,'node',Mobj.nVerts);
three_dimid=netcdf.defDim(nc,'three',3);
time_dimid=netcdf.defDim(nc,'time',netcdf.getConstant('NC_UNLIMITED'));
datestrlen_dimid=netcdf.defDim(nc,'DateStrLen',26);

% Space variables
lon_varid=netcdf.defVar(nc,'lon','NC_FLOAT',node_dimid);
netcdf.putAtt(nc,lon_varid,'long_name','nodal longitude');
netcdf.putAtt(nc,lon_varid,'units','degrees_easdt');

lat_varid=netcdf.defVar(nc,'lat','NC_FLOAT',node_dimid);
netcdf.putAtt(nc,lat_varid,'long_name','nodal latitude');
netcdf.putAtt(nc,lat_varid,'units','degrees_north');

lonc_varid=netcdf.defVar(nc,'lonc','NC_FLOAT',nele_dimid);
netcdf.putAtt(nc,lonc_varid,'long_name','zonal longitude');
netcdf.putAtt(nc,lonc_varid,'units','meters');

latc_varid=netcdf.defVar(nc,'latc','NC_FLOAT',nele_dimid);
netcdf.putAtt(nc,latc_varid,'long_name','zonal latitude');
netcdf.putAtt(nc,latc_varid,'units','meters');

nv_varid=netcdf.defVar(nc,'nv','NC_INT',[nele_dimid, three_dimid]);
netcdf.putAtt(nc,nv_varid,'long_name','nodes surrounding element');

% Time variables
time_varid=netcdf.defVar(nc,'time','NC_FLOAT',time_dimid);
netcdf.putAtt(nc,time_varid,'long_name','time');
netcdf.putAtt(nc,time_varid,'units','days since 1858-11-17 00:00:00');
netcdf.putAtt(nc,time_varid,'delta_t','0000-00-00 01:00:00')

netcdf.putAtt(nc,time_varid,'format','modified julian day (MJD)');
netcdf.putAtt(nc,time_varid,'time_zone','UTC');

times_varid=netcdf.defVar(nc,'Times','NC_CHAR',[datestrlen_dimid,time_dimid]);
netcdf.putAtt(nc,times_varid,'long_name','Calendar Date');
netcdf.putAtt(nc,times_varid,'format','String: Calendar Time');
netcdf.putAtt(nc,times_varid,'time_zone','UTC');

sst_varid = netcdf.defVar(nc, 'sst', 'NC_FLOAT', [node_dimid, time_dimid]);
netcdf.putAtt(nc, sst_varid, 'long_name', 'sea surface Temperature');
netcdf.putAtt(nc, sst_varid, 'units', 'Celsius Degree');
netcdf.putAtt(nc, sst_varid, 'grid', 'fvcom_grid');
netcdf.putAtt(nc, sst_varid, 'coordinates', Mobj.nativeCoords);
netcdf.putAtt(nc, sst_varid, 'type', 'data');

% End definitions
netcdf.endDef(nc);

% Put the easy ones in first.
netcdf.putVar(nc, nv_varid, Mobj.tri);
netcdf.putVar(nc,lon_varid,Mobj.lon);
netcdf.putVar(nc,lat_varid,Mobj.lat);
netcdf.putVar(nc,lonc_varid,Mobj.lonc);
netcdf.putVar(nc,latc_varid,Mobj.latc);
netcdf.putVar(nc,time_varid,0,ntimes,MJDtime);

nStringOut = char();
[nYr, nMon, nDay, nHour, nMin, nSec] = mjulian2greg(MJDtime);
for i=1:ntimes
    nDate = [nYr(i), nMon(i), nDay(i), nHour(i), nMin(i), nSec(i)];
    nStringOut = [nStringOut, sprintf('%04i/%02i/%02i %02i:%02i:%09.6f', nDate)];
end
netcdf.putVar(nc,times_varid,[0, 0], [26, ntimes], nStringOut);

netcdf.putVar(nc, sst_varid, [0, 0], [Mobj.nVerts, ntimes], sst)
fprintf('done.\n')

% Close the netCDF file(s)
netcdf.close(nc);

Mobj.assim.sst.data = sst;
Mobj.assim.sst.time = time;

% Plot sst if needed.
% for aa = 1:size(sst,2)
%     plot_field(Mobj,sst(:,aa))
%     pause
% end

if ftbverbose
    fprintf('end   : %s \n', subname)
end
