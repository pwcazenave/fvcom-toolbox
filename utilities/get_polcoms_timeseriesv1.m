function [polcoms]=get_polcoms_timeseriesv1(rootfname,ts_controlfile,Mobj,inputConf,tseries_dir,mm,polcoms)
%ts_controlfile='/users/modellers/rito/Models/MEDINA/tseries.MEDI29'
%tseries_dir='/data/perseus1/to_archive/suka_VECTORS/MEDI29_RA/physseries_MEDI29_2000';
%tseries_dir='/users/modellers/rito/Models/MEDINA/polcoms/physseries_MEDI29.2006';

% Run jobs on multiple workers if we have that functionality. Not sure if
% it's necessary, but check we have the Parallel Toolbox first.


%%
% Read Tseries file with number and locations of timeseries
nvars2 = 9;
var_names2 = {'u','v','temp','sal','aa','ak','qsq','al','iop'};
%%

fid = fopen(ts_controlfile);
nstations = textscan(fid, '%u',1);
C=  textscan(fid, '%u%u%f%f',nstations{1});
fclose(fid);
[lat,lon]=deal(C{3},C{4});
[polcoms_i,polcoms_j]=deal(C{1},C{2});
clear C
% convert lat and lon to utm
% Convert the small subdomain into cartesian coordinates.
tmpZone = regexpi(inputConf.utmZone,'\ ','split');
[tseries.x, tseries.y] = wgs2utm(lat(:), lon(:), str2double(char(tmpZone{1}(1))), 'N');

% Select points in the FVCOM domain (near the boundary as determined in
% Mobj
distance = abs(complex(tseries.x,tseries.y)-complex(nanmean(Mobj.x),nanmean(Mobj.y)));
dist_lim = mode(distance);
igood=find(distance < dist_lim*5);
%%

% build timeseries filenames for the time range under consideration
% read timeseries files. Make sure the correct variables in the files are
% read. We need a map of variables on the file
inputConf.zetUVfile=fullfile(tseries_dir,['zet_UBVB.',rootfname,'.',num2str(inputConf.modelYear),'.',num2str(mm,'%02d')]);
inputConf.PolcomsPoints=[polcoms_i(igood),polcoms_j(igood)];
% timeseries doesn't have information about the depth levels. Actual depths
% need to be calculated from total depth and scoord distribution.


% extract positions lat and lon at interest points
polcoms.bcidx=sub2ind(size(polcoms.latb),inputConf.PolcomsPoints(:,1),inputConf.PolcomsPoints(:,2));
latb=polcoms.latb(polcoms.bcidx);
lonb=polcoms.lonb(polcoms.bcidx);
[polcoms.bcxb, polcoms.bcyb] = wgs2utm(latb(:), lonb(:), str2double(char(tmpZone{1}(1))), 'N');
latu=polcoms.latu(polcoms.bcidx);
lonu=polcoms.lonu(polcoms.bcidx);
[polcoms.bcxu, polcoms.bcyu] = wgs2utm(latu(:), lonu(:), str2double(char(tmpZone{1}(1))), 'N');

% obtain depth levels for each station
% depth levels at each station need reconstructing because polcoms timeseries files do not
% include information on the depth levels.
[polcoms.xb, polcoms.yb] = wgs2utm(polcoms.latb(:), polcoms.lonb(:), str2double(char(tmpZone{1}(1))), 'N');
[polcoms.xu, polcoms.yu] = wgs2utm(polcoms.latu(:), polcoms.lonu(:), str2double(char(tmpZone{1}(1))), 'N');

fdb = TriScatteredInterp(polcoms.xb(:), polcoms.yb(:), polcoms.bathy(:), 'natural');
% interpolate bathymetry onto boundary points (b and u points)
polcoms.bchb=fdb(polcoms.bcxb,polcoms.bcyb);
polcoms.bchu=fdb(polcoms.bcxu, polcoms.bcyu);
% and onto fvcom bc positions (I don't think I need this)
% polcoms.hb=fdb(Mobj.x(oNodes),Mobj.y(oNodes));
% polcoms.hu=fdb(Mobj.xc(oElems),Mobj.yc(oElems));
polcoms.igood=igood;

%%


for ff=1:length(igood)
    fname =fullfile(tseries_dir,['physseries.',num2str(igood(ff)),'.',rootfname,'.',num2str(inputConf.modelYear),'.',num2str(mm,'%02d')]);
    cleanfile =fullfile(tseries_dir,'cleanfile');
    clean_statement=['sed ''s/^\**/0/g'' ' , fname,' > ',cleanfile];
    system(clean_statement);
    data = load(cleanfile);
    jday = data(:,1);data(:,1)=[];
    jday = reshape(jday,nvars2,[]);
    [~,ntimes]=size(jday);
    [~,ndepths]=size(data);
    
     jday = jday(1,:);jday = repmat(jday,[ndepths 1]);
    for nn=1:length(var_names2)
        polcoms.(var_names2{nn})(ff,:,:) = data(nn:length(var_names2):end,:)'./1000;
    end
     polcoms.jday(ff,:,:) = jday/24;
    
end
% Generate timerecord from filename and length of data
polcoms.time=datenum(inputConf.modelYear,mm,1):1/24:(size(polcoms.jday,3)-1)/24+datenum(inputConf.modelYear,mm,1);
inputConf.PolcomsLevs=ndepths;

% read zetUBVB file to extract surface elevation
[dumpstruct]=readzetUBVB(inputConf,ntimes);
polcoms=catstruct(dumpstruct,polcoms);
polcoms.ndepths=ndepths;
polcoms.ntimes=ntimes;
return
% make surface timeseries of 2d variables (zet, ub and vb)
% build timeseries matrix for each variable


% output as netcdf file in expected format for latest FVCOM 3.2

% Close the MATLAB pool if we opened it.

