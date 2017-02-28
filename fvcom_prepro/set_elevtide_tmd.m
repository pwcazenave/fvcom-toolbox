% Generate the required fields within Mobj for time series elevation
% forcing at open boundaries using the TMD (TPXO) functions.
%
% Note that TMD requires moving to the directory of the TMD toolbox. But
% this sould be OK.
%
% Rory O'Hara Murray, 2013-07-01
%
function [Mobj TMD_ConList] = set_elevtide_tmd(Mobj, dates_MJD)

%MyTitle = 'Julian FVCOM time series for open boundary from TPXO model using TMD';

% make sure an open boudary is defined
if(Mobj.nObs==0)
	warning('cannot setup spectral open boundary, there is no open boundary in the mesh struct')
	return
end;

% make sure lon and lat are defined
if(not(Mobj.have_lonlat))
    warning('cannot setup spectral open boundary, longitude and latitude are not defined')
    return
end

%for ob=1:Mobj.nObs % loop through each open boundary
ob = 1; %assume only one open bounary for the moment

BNid = Mobj.obc_nodes(ob,Mobj.obc_nodes(ob,:)>0); % get the ids of the nodes on the boundary
lat = Mobj.lat(BNid);
lon = Mobj.lon(BNid);

% generate 10 minute (600s) time series from TMD (TPXO)
%dates = [datenum(1993,09,1) datenum(1993,12,14)];
dates = dates_MJD([1 end]) + datenum('1858-11-17 00:00:00');
time = dates(1):1/24/6:dates(2);
time_MJD = time - datenum('1858-11-17 00:00:00');
model_file = 'DATA\Model_ES2008';
%tmp = which('tmd_tide_pred_2.m');
%a = strfind(tmp, 'tmd_tide_pred_2.m');
%cd(tmp(1:a-1))
% current_dir = pwd;
% cd([getenv('Hydro') '\Software\Matlab\ToolboxesExternal\TMD2.03\'])
[eta, TMD_ConList] = tmd_tide_pred_2(model_file, time, lat, lon, 'z');
% cd(current_dir);

%%
figure('position', [360   502   879   420])
%plot(NCOF_time, NCOF_eta(:,10), '-o', OTPS_time, OTPS_eta(:,10), OTPS_time, NCOF_eta2(:,10), '-')
t0 = time(1);
plot(time-t0, eta)
legend(gca, 'TMD (OTIS)', 4)
set(gca, 'yaxislocation', 'right')
xlabel('Days')
ylabel('Elevation (m)')
%title('September 2011')

%% save file in FVCOM 3 format (netCDF)
%write_FVCOM_julian(BNid,time,eta,filename_out_netCDF,MyTitle, 'timeformat', 'SDN') 

%% save time series to the Mobj structure
Mobj.surfaceElevation = eta';       % Makes surfaceElevation array with size (boundary nodes, time sereis)
Mobj.el_time = time_MJD;

return


