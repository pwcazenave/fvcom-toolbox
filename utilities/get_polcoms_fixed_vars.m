function [polcoms]=get_polcoms_fixed_vars(fileU,fileB,fileparams,filebathy,filescoord,ipexfile)
%
% fileU = '/data/perseus1/to_archive/suka_VECTORS/MEDI29_RA/daily/dmeanUVT.MEDI29.RA.2000.01.U.nc'
% fileB =  '/data/perseus1/to_archive/suka_VECTORS/MEDI29_RA/daily/dmeanUVT.MEDI29.RA.2000.01.B.nc'
% filebio='/data/perseus1/to_archive/suka_VECTORS/MEDI29_RA/daily/dmeaneco3d.MEDI29.RA.2000.01.nc'
% ipexfile='/users/modellers/rito/Models/MEDINA/polcoms/zet_UBVB.MEDI29.RA.2009.12'
% filescoord='/users/modellers/rito/Models/MEDINA/polcoms/scoord_params.dat'
% fileparams='/users/modellers/rito/Models/MEDINA/polcoms/MEDI29.parameters'
% filebathy='/users/modellers/rito/Models/MEDINA/polcoms/MEDI29.bathy'
% Read Tseries file with number and locations of timeseries
varlistU = {'lon', 'lat'};
varlistB = {'lon', 'lat'};
% varlistbio = {'depth', 'pdepth'};

% read s_coord file from polcoms setup


polcoms.scoord=read_scoord_params(filescoord);
% Read domain parameters
polcoms.params=read_polcoms_params(fileparams);
% read ipexu and ipexb
%fid is for the input file
fid=fopen(ipexfile,'r','n');
dump= fread(fid,1,'int32');
polcoms.iesub= fread(fid,1,'int32');
polcoms.jesub= fread(fid,1,'int32');
polcoms.n= fread(fid,1,'int32');
polcoms.npsea= fread(fid,1,'int32');
dump = fread(fid,2,'int32');
polcoms.isea= fread(fid,polcoms.npsea,'int32');
dump = fread(fid,2,'int32');
polcoms.jsea= fread(fid,polcoms.npsea,'int32');
dump = fread(fid,2,'int32');

l= fread(fid,1,'int32');
m= fread(fid,1,'int32');
n= fread(fid,1,'int32');
polcoms.npusea= fread(fid,1,'int32');
dump = fread(fid,2,'int32');
polcoms.iusea= fread(fid,polcoms.npusea,'int32');
dump = fread(fid,2,'int32');
polcoms.jusea= fread(fid,polcoms.npusea,'int32');
dump = fread(fid,2,'int32');
fclose(fid);
% read bathymetry and calculate scoords.
polcoms.bathy = textread(filebathy,'%f','delimiter',' ');
polcoms.bathy = reshape(polcoms.bathy,polcoms.iesub,polcoms.jesub);
% 
pcU = get_POLCOMS_netCDF(fileU, varlistU);
pcB = get_POLCOMS_netCDF(fileB, varlistB);
% pcbio=get_POLCOMS_netCDF(filebio, varlistbio);
[polcoms.lonb,polcoms.latb]=meshgrid(pcB.lon.data,pcB.lat.data);
[polcoms.lonu,polcoms.latu]=meshgrid(pcU.lon.data,pcU.lat.data);

polcoms.lonb=polcoms.lonb';
polcoms.latb=polcoms.latb';
polcoms.lonu=polcoms.lonu';
polcoms.latu=polcoms.latu';
[polcoms]=calc_scoord(polcoms);
return
