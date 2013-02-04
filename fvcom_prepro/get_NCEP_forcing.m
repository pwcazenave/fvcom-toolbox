function data = get_NCEP_forcing(Mobj, modelTime)
% Get the required parameters from NCEP OPeNDAP data to force FVCOM
% (through any of Casename_wnd.nc, Casename_sst.nc, Casename_hfx.nc
% or Casename_pre_evap.nc).
%
% data = get_NCEP_forcing(Mobj, modelTime)
%
% DESCRIPTION:
%   Using OPeNDAP, extract the necessary parameters to create an FVCOM
%   forcing file. Requires the air_sea toolbox and the OPeNDAP toolbox (see
%   below for where to get them).
%
% INPUT:
%   Mobj - MATLAB mesh object
%   modelTime - Modified Julian Date start and end times
%
% OUTPUT:
%   data - struct of the data necessary to force FVCOM. These can be
%   interpolated onto an unstructured grid in Mobj using
%   grid2fvcom_U10V10.m.
%
% The parameters which can be obtained from the NCEP data are:
%     - u wind component (uwnd)
%     - v wind component (vwnd)
%     - Net longwave radiation surface (nlwrs)
%     - Net shortwave radiation surface (nswrs)
%     - Air temperature (air)
%     - Relative humidity (rhum)
%     - Precipitation rate (prate)
%     - Sea level pressure (slp)
%     - Latent heat flux (lhtfl)
%     - Surface heat flux (shtfl)
%     - Potential evaporation rate (pevpr)
%
% In addition to these, the momentum flux is calculated from wind data.
% Precipitation is converted from kg/m^2/s to m/s. Evaporation is
% calculated from the mean daily latent heat net flux (lhtfl) at the
% surface.
%
% REQUIRES:
%   The air_sea toolbox:
%       http://woodshole.er.usgs.gov/operations/sea-mat/air_sea-html/index.html
%   The OPeNDAP toolbox:
%       http://www.opendap.org/pub/contributed/source/ml-toolbox/
%
%
% Author(s)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%   Ricardo Torres (Plymouth Marine Laboratory)
%
% Revision history:
%   2012-10-31 First version based on get_NCEP_L4.m.
%
%==========================================================================

subname = 'get_NCEP_forcing';

global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end

% Get the extent of the model domain (in spherical)
if ~Mobj.have_lonlat
    error('Need spherical coordinates to extract the forcing data')
else
    % Add a buffer of one grid cell in latitude and two in longitude to
    % make sure the model domain is fully covered by the extracted data.
    [dx, dy] = deal(2.5, 2.5); % NCEP resolution in degrees
    extents = [min(Mobj.lon(:))-(2*dx), max(Mobj.lon(:))+(2*dx), min(Mobj.lat(:))-dy, max(Mobj.lat(:))+dy];
end

if modelTime(end) - modelTime(1) > 365
    error('Can''t (yet) process more than a year at a time.')
end

yearStart = mjulian2greg(modelTime(1));
yearEnd = mjulian2greg(modelTime(end));
if yearEnd ~= yearStart
    error('Can''t (yet) process across a year boundary.')
else
    year = yearEnd;
end

% Set up a struct of the NCEP remote locations in which we're interested.
ncep.uwnd   = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface_gauss/uwnd.10m.gauss.',num2str(year),'.nc'];
ncep.vwnd   = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface_gauss/vwnd.10m.gauss.',num2str(year),'.nc'];
ncep.nlwrs  = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface_gauss/nlwrs.sfc.gauss.',num2str(year),'.nc'];
ncep.nswrs  = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface_gauss/nswrs.sfc.gauss.',num2str(year),'.nc'];
ncep.air    = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface_gauss/air.2m.gauss.',num2str(year),'.nc'];
ncep.rhum   = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface/rhum.sig995.',num2str(year),'.nc'];
ncep.prate  = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface_gauss/prate.sfc.gauss.',num2str(year),'.nc'];
ncep.slp    = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface/slp.',num2str(year),'.nc'];
ncep.lhtfl  = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface_gauss/lhtfl.sfc.gauss.',num2str(year),'.nc'];
ncep.shtfl  = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface_gauss/shtfl.sfc.gauss.',num2str(year),'.nc'];
ncep.pevpr  = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface_gauss/pevpr.sfc.gauss.',num2str(year),'.nc'];

fields = fieldnames(ncep);

for aa=1:length(fields)
    data.(fields{aa}).data = [];
    data.(fields{aa}).time = [];
    data.(fields{aa}).lat = [];
    data.(fields{aa}).lon = [];
    data_attributes.(fields{aa}) = [];

    % Get the data time and convert to Modified Julian Day.
    data_time = loaddap([ncep.(fields{aa}),'?time']);
    data_attributes.(fields{aa}) = loaddap('-A',[ncep.(fields{aa})]);
    timevec = datevec((data_time.time)/24+365);
    data.time = greg2mjulian(timevec(:,1), timevec(:,2), timevec(:,3), ...
        timevec(:,4), timevec(:,5), timevec(:,6));
    % Clip the time to the given range
    data_time_mask = data.time >= modelTime(1) & data.time <= modelTime(end);
    data_time_idx = 1:size(data.time,1);
    data_time_idx = data_time_idx(data_time_mask);
    data.time = data.time(data_time_mask);

    % Check the times
    %[yyyy,mm,dd,hh,MM,ss] = mjulian2greg(data.time(1))
    %[yyyy,mm,dd,hh,MM,ss] = mjulian2greg(data.time(end))

    % Clip the data to the model domain
    data_lon = loaddap([ncep.(fields{aa}),'?lon']);
    % If the extents are negative in longitude, we need to extract the NCEP
    % data in two goes, once for the end of the grid (west of Greenwich),
    % once for the beginning (east of Greenwich), and then stick the two
    % bits together.
    clear index_lon index_lat
    if extents(1) < 0 && extents(2) < 0
        % This is OK, we can just shunt the values by 360.
        extents(1) = extents(1) + 360;
        extents(2) = extents(2) + 360;
        index_lon = find(data_lon.lon > extents(1) & data_lon.lon < extents(2));
    elseif extents(1) < 0 && extents(2) > 0
        % This is the tricky one. We'll do two passes to extract the
        % western chunk first (extents(1)+360 to 360), then the eastern
        % chunk (0-extents(2))
        index_lon{1} = find(data_lon.lon >= extents(1) + 360);
        index_lon{2} = find(data_lon.lon <= extents(2));
    else
        % Dead easy, we're in the eastern hemisphere, so nothing too
        % strenuous here
        index_lon = find(data_lon.lon > extents(1) & data_lon.lon < extents(2));
    end

    % Latitude is much more straightforward
    data_lat = loaddap([ncep.(fields{aa}),'?lat']);
    index_lat = find(data_lat.lat > extents(3) & data_lat.lat < extents(4));

    % Get the data
    if iscell(index_lon)
        % We need to do each half and merge them
        eval(['data1_west.(fields{aa}) = loaddap(''', ncep.(fields{aa}),'?',...
            fields{aa},'[', num2str(min(data_time_idx)-1),':',...
            num2str(max(data_time_idx)-1), '][',...
            num2str(min(index_lat)-1), ':', num2str(max(index_lat)-1),...
            '][', num2str(min(index_lon{1})-1), ':',...
            num2str(length(data_lon.lon)-1), ']'');']);
        eval(['data1_east.(fields{aa}) = loaddap(''', ncep.(fields{aa}),'?',...
            fields{aa}, '[', num2str(min(data_time_idx)-1),':',...
            num2str(max(data_time_idx)-1), '][',...
            num2str(min(index_lat)-1), ':', num2str(max(index_lat)-1),...
            '][', '0', ':', num2str(max(index_lon{2})-1), ']'');']);
        % Merge the two sets of data together
        structfields = fieldnames(data1_west.(fields{aa}).(fields{aa}));
        for ii=1:length(structfields)
            switch structfields{ii}
                case 'lon'
                    % Only the longitude and the actual data need sticking
                    % together, but each must be done along a different
                    % axis (lon is a vector, the data is an array).
                    data1.(fields{aa}).(fields{aa}).(structfields{ii}) = ...
                        [data1_west.(fields{aa}).(fields{aa}).(structfields{ii});data1_east.(fields{aa}).(fields{aa}).(structfields{ii})];
                case fields{aa}
                    % This is the actual data
                    data1.(fields{aa}).(fields{aa}).(structfields{ii}) = ...
                        [data1_west.(fields{aa}).(fields{aa}).(structfields{ii}),data1_east.(fields{aa}).(fields{aa}).(structfields{ii})];
                otherwise
                    % Assume the data are the same in both arrays. A simple
                    % check of the range of values in the difference
                    % between the two arrays should show whether they're
                    % the same or not. If they are, use the western values,
                    % otherwise, warn about the differences. It might be
                    % the data are relatively unimportant anyway (i.e. not
                    % used later on).
                    try
                        tdata = data1_west.(fields{aa}).(fields{aa}).(structfields{ii}) - data1_east.(fields{aa}).(fields{aa}).(structfields{ii});
                        if range(tdata(:)) == 0
                            % They're the same data
                            data1.(fields{aa}).(fields{aa}).(structfields{ii}) = ...
                                data1_west.(fields{aa}).(fields{aa}).(structfields{ii});
                        else
                            warning('Unexpected data field and the west and east halves don''t match. Skipping.')
                        end
                    catch
                        warning('Unexpected data field and the west and east halves don''t match. Skipping.')
                    end
                    clear tdata
            end
        end
    else
        % We have a straightforward data extraction
        eval(['data1.(fields{aa}) = loaddap(''', ncep.(fields{aa}),'?',...
            fields{aa}, '[', num2str(min(data_time_idx)-1),':',...
            num2str(max(data_time_idx)-1), '][',...
            num2str(min(index_lat)-1), ':', num2str(max(index_lat)-1),...
            '][', num2str(min(index_lon)-1), ':',...
            num2str(max(index_lon)-1), ']'');']);
    end

    datatmp = squeeze(data1.(fields{aa}).(fields{aa}).(fields{aa}));
    datatmp = (datatmp * data_attributes.(fields{aa}).(fields{aa}).scale_factor) + data_attributes.(fields{aa}).(fields{aa}).add_offset;

    data.(fields{aa}).data = cat(1, data.(fields{aa}).data, datatmp);
    data.(fields{aa}).time = cat(1, data.(fields{aa}).time, squeeze(data1.(fields{aa}).(fields{aa}).time));
    data.(fields{aa}).lat = squeeze(data1.(fields{aa}).(fields{aa}).lat);
    data.(fields{aa}).lon = squeeze(data1.(fields{aa}).(fields{aa}).lon);
end

% Now we have some data, we need to create some additional parameters
% required by FVCOM.

% Convert precipitation from kg/m^2/s to m/s (required by FVCOM) by
% dividing by freshwater density (kg/m^3).
data.prate.data = data.prate.data/1000;

% Evaporation can be approximated by:
%
%   E(m/s) = lhtfl/Llv/rho
%
% where:
%
%   lhtfl   = "Mean daily latent heat net flux at the surface"
%   Llv     = Latent heat of vaporization (approx to 2.5*10^6 J kg^-1)
%   rho     = 1025 kg/m^3
%
Llv = 2.5*10^6;
rho = 1025; % using a typical value for seawater.
Et = data.lhtfl.data/Llv/rho;
data.P_E.data = data.prate.data-Et;

% Calculate the momentum flux
WW = data.uwnd.data + data.vwnd.data * 1i;
data.tau.data = stresslp(abs(WW),10);
[data.tx.data,data.ty.data] = wstress(data.uwnd.data,data.vwnd.data,10);
data.tx.data=reshape(data.tx.data*0.1, size(data.uwnd.data)); % dyn/cm^2 to N/m^2
data.ty.data=reshape(data.ty.data*0.1, size(data.uwnd.data)); % dyn/cm^2 to N/m^2

% Get the fields we need for the subsequent interpolation
data.lon = data.uwnd.lon;
data.lon(data.lon > 180) = data.lon(data.lon > 180) - 360;
data.lat = data.uwnd.lat;

% Have a look at some data.
% [X, Y] = meshgrid(data.lon, data.lat);
% for i=1:size(data.uwnd.data, 3)
%     figure(1)
%     clf
%     uv = sqrt(data.uwnd.data(:, :, i).^2 + data.vwnd.data(:, :, i).^2);
%     pcolor(X, Y, uv)
%     shading flat
%     axis('equal','tight')
%     pause(0.1)
% end
