function write_nesting_bdy_file(Mobj, MJD, OutFile, MyTitle)
% Write an FVCOM time series forcing file
%
% write_wave_bdy_file(Mobj, MJD, OutFile, MyTitle)
%
% DESCRIPTION:
%    Write an FVCOM NetCDF  nesting boundary file
%
% INPUT:
%   Mobj = Matlab mesh object.
%   MJD = list of modified Modified Julian Dates of size [times] (defined
%         as unlimited in the NetCDF file).
%   OutFile = name of NetCDF file.
%   MyTitle = casename title, written as global attribute of NetCDF file.
%
%
% OUTPUT:
%    WaveFile, A NetCDF FVCOM-SWAVE wave forcing file
%
% EXAMPLE USAGE
%    write_FVCOM_elevtide(Mobj, MJD, '/tmp/waves.nc, 'Shelf tides')
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Karen Thurston (National Oceanography Centre Liverpool)
%    DarrenPrice (based upon above authors)
% Revision history
%    2012-08-08 (PWC) First version.
%    2012-11-14 (PWC) Updated to expect Modified Julian Day rather than
%    doing the conversion in here. Also put the pieces in set_elevtide in
%    here to simplify the process of writing out an elevation input file.
%    2012-12-04 (KJT) Updated to use surface elevation and open boundary
%    nodes from Mobj
%    2013-11-05 (DMP) Updated to use wave conditions and open boundary
%    2014-01-20 (HKJ) Corrected for wave conditions at open boundary
%    2014-06-03 (DMP) Updated to use Z, U,V,T,S conditions and open boundary
%
%==========================================================================

global ftbverbose
report = false;
if(ftbverbose); report = true; end
subname = 'write_nesting_bdy_file';
if(report); fprintf('\n'); end
if(report); fprintf(['begin : ' subname '\n']); end

%{
% ---  extract from ... write_elevation..
% Get a list of the open boundary nodes. Transpose Mobj.obc_nodes so the
% order of the boundary nodes is preserved.
tmpObcNodes = Mobj.nnodesID';
% Flip it back so it's the same shape as it would have been using the old
% code.
ObcNodes = tmpObcNodes(tmpObcNodes~=0)';
%}

ObcNodes = Mobj.nnodesID;
nElems  = Mobj.nElems;

%--------------------------------------------------------------------------
% Sanity check on input and dimensions
%--------------------------------------------------------------------------
nTimes = numel(MJD);
if(report); fprintf('Number of time steps %d\n',nTimes); end

nObcs = numel(ObcNodes);
if(report); fprintf('Number of Open Boundary Nodes %d\n',nObcs); end

[chk1, chk2] = size(Mobj.surfaceElevation);
if nObcs ~= chk1 || nTimes ~= chk2
    fprintf('surface Elevation dimensions do not match time series and number of boundary nodes.\n')
    fprintf('surface Elevation nodes and time sizes: (%d, %d)\n', chk1, chk2)
    fprintf('Boundary nodes size: %d\n', nObcs)
    fprintf('Times size: %d\n', nTimes)
	error('Input data sizes do not match. Check and try again.');
end

%--------------------------------------------------------------------------
% Dump the file
%--------------------------------------------------------------------------

nc=netcdf.create(OutFile,'clobber');

% define global attributes
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title',MyTitle)
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type','FVCOM TIME SERIES NESTING FORCING FILE')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history', sprintf('File created using %s from the MATLAB fvcom-toolbox', subname))
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.0')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'CoordinateSystem','Spherical')
%netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'CoordinateSystem','none')
%disp('need to sort out coordinate system in write_nesting_bdy_file.m');




% define dimensions
%hkj nobc_dimid=netcdf.defDim(nc,'nobc',nObcs);
nlevels = Mobj.nsiglev;
nlayers = nlevels - 1;
nele_dimid  = netcdf.defDim(nc,'nele',nElems);
nobc_dimid  = netcdf.defDim(nc,'node',nObcs);
siglay_dimid= netcdf.defDim(nc,'siglay',nlayers);
siglev_dimid= netcdf.defDim(nc,'siglev',nlevels);
three_dimid = netcdf.defDim(nc,'three',3);
time_dimid  = netcdf.defDim(nc,'time',netcdf.getConstant('NC_UNLIMITED'));
one_dimid = netcdf.defDim(nc,'one',1);
date_str_len_dimid=netcdf.defDim(nc,'DateStrLen',26);



% define variables and attributes
disp('writing nprocs attributes');
nprocs_varid=netcdf.defVar(nc,'nprocs','NC_INT',one_dimid);
netcdf.putAtt(nc,nprocs_varid,'long_name','number of processors');

disp('writing npartition attributes');
npart_varid=netcdf.defVar(nc,'partition','NC_INT',nele_dimid);
netcdf.putAtt(nc,npart_varid,'long_name','partition');

disp('writing x attributes');
x_varid=netcdf.defVar(nc,'x','NC_FLOAT',nobc_dimid);
netcdf.putAtt(nc,x_varid,'long_name','nodal x-coordinate');
netcdf.putAtt(nc,x_varid,'units','meters');

disp('writing y attributes');
y_varid=netcdf.defVar(nc,'y','NC_FLOAT',nobc_dimid);
netcdf.putAtt(nc,y_varid,'long_name','nodal y-coordinate');
netcdf.putAtt(nc,y_varid,'units','meters');

disp('writing lon attributes');
lon_varid=netcdf.defVar(nc,'lon','NC_FLOAT',nobc_dimid);
netcdf.putAtt(nc,lon_varid,'long_name','nodal longitude');
netcdf.putAtt(nc,lon_varid,'standard_name','longitude');
netcdf.putAtt(nc,lon_varid,'units','degrees_east');

disp('writing lat attributes');
lat_varid=netcdf.defVar(nc,'lat','NC_FLOAT',nobc_dimid);
netcdf.putAtt(nc,lat_varid,'long_name','nodal latitude');
netcdf.putAtt(nc,lat_varid,'standard_name','latitude');
netcdf.putAtt(nc,lat_varid,'units','degrees_north');

disp('writing xc attributes');
xc_varid=netcdf.defVar(nc,'xc','NC_FLOAT',nele_dimid);
netcdf.putAtt(nc,xc_varid,'long_name','zonal x-coordinate');
netcdf.putAtt(nc,xc_varid,'units','meters');

disp('writing yc attributes');
yc_varid=netcdf.defVar(nc,'yc','NC_FLOAT',nele_dimid);
netcdf.putAtt(nc,yc_varid,'long_name','zonal y-coordinate');
netcdf.putAtt(nc,yc_varid,'units','meters');

disp('writing lonc attributes');
lonc_varid=netcdf.defVar(nc,'lonc','NC_FLOAT',nele_dimid);
netcdf.putAtt(nc,lonc_varid,'long_name','zonal longitude');
netcdf.putAtt(nc,lonc_varid,'standard_name','longitude');
netcdf.putAtt(nc,lonc_varid,'units','degrees_east');

disp('writing latc attributes');
latc_varid=netcdf.defVar(nc,'latc','NC_FLOAT',nele_dimid);
netcdf.putAtt(nc,latc_varid,'long_name','zonal latitude');
netcdf.putAtt(nc,latc_varid,'standard_name','latitude');
netcdf.putAtt(nc,latc_varid,'units','degrees_north');

disp('writing siglay attributes');
siglay_varid=netcdf.defVar(nc,'siglay','NC_FLOAT',[nobc_dimid, siglay_dimid ]);
netcdf.putAtt(nc,siglay_varid,'long_name','Sigma Layers');
netcdf.putAtt(nc,siglay_varid,'standard_name','ocean_sigma/general_coordinate');
netcdf.putAtt(nc,siglay_varid,'positive','up');
netcdf.putAtt(nc,siglay_varid,'valid_min','-1.f');
netcdf.putAtt(nc,siglay_varid,'valid_max','0.f');
netcdf.putAtt(nc,siglay_varid,'formula_terms','sigma: siglay eta: zeta depth: h');

disp('writing siglev attributes');
siglev_varid=netcdf.defVar(nc,'siglev','NC_FLOAT',[nobc_dimid,siglev_dimid ]);
netcdf.putAtt(nc,siglev_varid,'long_name','Sigma Levels');
netcdf.putAtt(nc,siglev_varid,'standard_name','ocean_sigma/general_coordinate');
netcdf.putAtt(nc,siglev_varid,'positive','up');
netcdf.putAtt(nc,siglev_varid,'valid_min','-1.f');
netcdf.putAtt(nc,siglev_varid,'valid_max','0.f');
netcdf.putAtt(nc,siglev_varid,'formula_terms','sigma: siglay eta: zeta depth: h');

disp('writing h attributes');
h_varid=netcdf.defVar(nc,'h','NC_FLOAT',nobc_dimid);
netcdf.putAtt(nc,h_varid,'long_name','Bathymetry');
netcdf.putAtt(nc,h_varid,'standard_name','sea_floor_depth_below_geoid');
netcdf.putAtt(nc,h_varid,'units','m');
netcdf.putAtt(nc,h_varid,'positive','down');
netcdf.putAtt(nc,h_varid,'grid','Bathymetry_Mesh');
netcdf.putAtt(nc,h_varid,'coordinates','lat lon');
netcdf.putAtt(nc,h_varid,'type','data');

disp('writing nv attributes');
% nv_varid=netcdf.defVar(nc,'nv','NC_INT',[three_dimid, nele_dimid]);
nv_varid=netcdf.defVar(nc,'nv','NC_INT',[nele_dimid, three_dimid]);
netcdf.putAtt(nc,nv_varid,'long_name','nodes surrounding element');

%hkj test - unsure if this is required
disp('writing obc_nodes attributes');
nobc_varid=netcdf.defVar(nc,'obc_nodes','NC_INT',nobc_dimid);
netcdf.putAtt(nc,nobc_varid,'long_name','Open Boundary Node Number');
netcdf.putAtt(nc,nobc_varid,'grid','obc_grid');

disp('writing iint attributes');
iint_varid=netcdf.defVar(nc,'iint','NC_INT',time_dimid);
netcdf.putAtt(nc,iint_varid,'long_name','internal mode iteration number');

disp('writing time attributes');
time_varid=netcdf.defVar(nc,'time','NC_FLOAT',time_dimid);
netcdf.putAtt(nc,time_varid,'long_name','time');
netcdf.putAtt(nc,time_varid,'units','days since 1858-11-17 00:00:00');
netcdf.putAtt(nc,time_varid,'format','modified julian day (MJD)');
netcdf.putAtt(nc,time_varid,'time_zone','UTC');

disp('writing Itime attributes');
itime_varid=netcdf.defVar(nc,'Itime','NC_INT',time_dimid);
netcdf.putAtt(nc,itime_varid,'units','days since 1858-11-17 00:00:00');
netcdf.putAtt(nc,itime_varid,'format','modified julian day (MJD)');
netcdf.putAtt(nc,itime_varid,'time_zone','UTC');

disp('writing Itime2 attributes');
itime2_varid=netcdf.defVar(nc,'Itime2','NC_INT',time_dimid);
netcdf.putAtt(nc,itime2_varid,'units','msec since 00:00:00');
netcdf.putAtt(nc,itime2_varid,'time_zone','UTC');

disp('writing Times attributes');
Times_varid=netcdf.defVar(nc,'Times','NC_CHAR',[date_str_len_dimid, time_dimid]);
netcdf.putAtt(nc,Times_varid,'time_zone','UTC');

disp('writing salinity attributes');
salinity_varid=netcdf.defVar(nc,'salinity','NC_FLOAT',[nobc_dimid, siglay_dimid, time_dimid]);
netcdf.putAtt(nc,salinity_varid,'long_name','salinity');
netcdf.putAtt(nc,salinity_varid,'standard_name','sea_water_salinity');
netcdf.putAtt(nc,salinity_varid,'units','1e-3');
netcdf.putAtt(nc,salinity_varid,'grid','fvcom_grid');
netcdf.putAtt(nc,salinity_varid,'coordinates','lat lon');
netcdf.putAtt(nc,salinity_varid,'type','data');

disp('writing temperature attributes');
temperature_varid=netcdf.defVar(nc,'temp','NC_FLOAT',[nobc_dimid,siglay_dimid, time_dimid]);
netcdf.putAtt(nc,temperature_varid,'long_name','temperature');
netcdf.putAtt(nc,temperature_varid,'standard_name','sea_water_temperature');
netcdf.putAtt(nc,temperature_varid,'units','degrees_C');
netcdf.putAtt(nc,temperature_varid,'grid','fvcom_grid');
netcdf.putAtt(nc,temperature_varid,'coordinates','lat lon');
netcdf.putAtt(nc,temperature_varid,'type','data');

disp('writing u velocity attributes');
u_varid=netcdf.defVar(nc,'u','NC_FLOAT',[nele_dimid,siglay_dimid, time_dimid]);
netcdf.putAtt(nc,u_varid,'long_name','Eastward Water Velocity');
netcdf.putAtt(nc,u_varid,'units','meters s-1');
netcdf.putAtt(nc,u_varid,'grid','fvcom_grid');
netcdf.putAtt(nc,u_varid,'type','data');

disp('writing v velocity attributes');
v_varid=netcdf.defVar(nc,'v','NC_FLOAT',[nele_dimid,siglay_dimid, time_dimid]);
netcdf.putAtt(nc,v_varid,'long_name','Northward Water Velocity');
netcdf.putAtt(nc,v_varid,'units','meters s-1');
netcdf.putAtt(nc,v_varid,'grid','fvcom_grid');
netcdf.putAtt(nc,v_varid,'type','data');

disp('writing vertical velocity attributes');
omega_varid=netcdf.defVar(nc,'omega','NC_FLOAT',[nobc_dimid,siglev_dimid, time_dimid]);
netcdf.putAtt(nc,omega_varid,'long_name','hydro static vertical velocity');
netcdf.putAtt(nc,omega_varid,'units','meters s-1');
netcdf.putAtt(nc,omega_varid,'grid','fvcom_grid');
netcdf.putAtt(nc,omega_varid,'coordinates','lat lon');
netcdf.putAtt(nc,omega_varid,'type','data');

disp('writing vertical velocity attributes');
hyw_varid=netcdf.defVar(nc,'hyw','NC_FLOAT',[nobc_dimid,siglev_dimid, time_dimid]);
netcdf.putAtt(nc,hyw_varid,'long_name','hydro static vertical velocity');
netcdf.putAtt(nc,hyw_varid,'units','meters s-1');
netcdf.putAtt(nc,hyw_varid,'grid','fvcom_grid');
netcdf.putAtt(nc,hyw_varid,'coordinates','lat lon');
netcdf.putAtt(nc,hyw_varid,'type','data');

disp('writing depth averaged u velocity attributes');
ua_varid=netcdf.defVar(nc,'ua','NC_FLOAT',[nele_dimid, time_dimid]);
netcdf.putAtt(nc,ua_varid,'long_name','Vertically Averaged x-velocity');
netcdf.putAtt(nc,ua_varid,'units','meters s-1');
netcdf.putAtt(nc,ua_varid,'grid','fvcom_grid');
netcdf.putAtt(nc,ua_varid,'type','data');

disp('writing depth averaged v velocity attributes');
va_varid=netcdf.defVar(nc,'va','NC_FLOAT',[nele_dimid, time_dimid]);
netcdf.putAtt(nc,va_varid,'long_name','Vertically Averaged y-velocity');
netcdf.putAtt(nc,va_varid,'units','meters s-1');
netcdf.putAtt(nc,va_varid,'grid','fvcom_grid');
netcdf.putAtt(nc,va_varid,'type','data');

disp('writing water level attributes');
zeta_varid=netcdf.defVar(nc,'zeta','NC_FLOAT',[nobc_dimid, time_dimid]);
netcdf.putAtt(nc,zeta_varid,'long_name','Water Surface Elevation');
netcdf.putAtt(nc,zeta_varid,'units','meters');
netcdf.putAtt(nc,zeta_varid,'positive','up');
netcdf.putAtt(nc,zeta_varid,'standard_name','sea_surface_elevation');
netcdf.putAtt(nc,zeta_varid,'grid','SSH_Mesh');
netcdf.putAtt(nc,zeta_varid,'coordinates','lat lon');
netcdf.putAtt(nc,zeta_varid,'type','data');


% for type 3 nesting need to include weighting
if (Mobj.NEST_TYPE=='TYPE3')

    disp('writing node weighting attributes');
    node_weight_varid=netcdf.defVar(nc,'weight_node','NC_FLOAT',[nobc_dimid time_dimid]);
    netcdf.putAtt(nc,node_weight_varid,'long_name','nesting node weighting');
    netcdf.putAtt(nc,node_weight_varid,'units','none');
    netcdf.putAtt(nc,node_weight_varid,'grid','fvcom_grid');
    netcdf.putAtt(nc,node_weight_varid,'coordinates','lat lon');
    netcdf.putAtt(nc,node_weight_varid,'type','data');

        disp('writing element weighting attributes');
    elem_weight_varid=netcdf.defVar(nc,'weight_cell','NC_FLOAT',[nele_dimid time_dimid]);
    netcdf.putAtt(nc,elem_weight_varid,'long_name','nesting element weighting');
    netcdf.putAtt(nc,elem_weight_varid,'units','none');
    netcdf.putAtt(nc,elem_weight_varid,'grid','fvcom_grid');
    netcdf.putAtt(nc,node_weight_varid,'coordinates','latc lonc');
    netcdf.putAtt(nc,elem_weight_varid,'type','data');


end

% end definitions
netcdf.endDef(nc);

%------------------------------------------------------------
% write data
%------------------------------------------------------------
disp('writing nprocs - not used??');
nprocs(1)= 3;
netcdf.putVar(nc,nprocs_varid,nprocs);

disp('writing partition - not used??');
partition = ones(1,nElems);
netcdf.putVar(nc,npart_varid,partition);

disp('writing co-ordinate data at nodal points: x, y, lon, lat');
netcdf.putVar(nc,x_varid,Mobj.x)
netcdf.putVar(nc,y_varid,Mobj.y)
netcdf.putVar(nc,lon_varid,Mobj.lon)
netcdf.putVar(nc,lat_varid,Mobj.lat)

disp('writing co-ordinate data at element centres: xc, yc, lonc, latc');
netcdf.putVar(nc,xc_varid,Mobj.xc)
netcdf.putVar(nc,yc_varid,Mobj.yc)
netcdf.putVar(nc,lonc_varid,Mobj.lonc)
netcdf.putVar(nc,latc_varid,Mobj.latc)

disp('writing vertical schematization data: sigma layers, sigma levels');
inc = 1./real(nlayers);
siglev = 0:-inc:-1;
for i=1:nlevels
	Msiglev(i,1:nObcs) = siglev(i);
end;
for i=1:nlayers
	Msiglay(i,1:nObcs) = mean(siglev(i:i+1));
end;
%{
Msiglay = -1*ones(nlayers,nObcs);
tlayer = 1/nlayers;
temp   = 0.5*tlayer;
for i=1:nlayers
    Msiglay(i,:)= temp*Msiglay(i,:);
    temp        = temp + tlayer;
end

Msiglev = -1*ones(nlevels,nObcs);
temp    = 0.0;
for i=1:nlevels
    Msiglev(i,:)= temp*Msiglev(i,:);
    temp        = temp + tlayer;
end
%}
netcdf.putVar(nc,siglay_varid,Msiglay)
netcdf.putVar(nc,siglev_varid,Msiglev)

disp('writing h_varid data');
netcdf.putVar(nc,h_varid,Mobj.h);

disp('writing nv_varid data');
Mobj.nv = Mobj.nv';
netcdf.putVar(nc,nv_varid,Mobj.nv);

%hkj test - doubtful if this is required
disp('writing nobc_varid data');
netcdf.putVar(nc,nobc_varid,ObcNodes);

disp('writing iint_varid data');
netcdf.putVar(nc,iint_varid,0,nTimes,1:nTimes);

disp('writing time_varid data');
netcdf.putVar(nc,time_varid,0,nTimes,MJD);

disp('writing itime_varid data');
netcdf.putVar(nc,itime_varid,floor(MJD));

disp('writing itime2_varid data');
netcdf.putVar(nc,itime2_varid,0,nTimes,mod(MJD,1)*24*3600*1000);  %

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

disp('writing Times_varid data');
netcdf.putVar(nc,Times_varid,nStringOut);

% disp('writing Times_varid data');
% netcdf.putVar(nc,Times_varid,Mobj.Times);  %Mobj.Times);

disp('writing salinity_varid data');
netcdf.putVar(nc,salinity_varid,Mobj.salinity);

disp('writing temperature_varid data');
netcdf.putVar(nc,temperature_varid,Mobj.temperature);

disp('writing u_varid data');
netcdf.putVar(nc,u_varid,Mobj.u);

disp('writing v_varid data');
netcdf.putVar(nc,v_varid,Mobj.v);

disp('writing ua_varid data');
netcdf.putVar(nc,ua_varid,Mobj.daUvel);

disp('writing va_varid data');
netcdf.putVar(nc,va_varid,Mobj.daVvel);

disp('writing omega_varid data');
netcdf.putVar(nc,omega_varid,Mobj.w);

disp('writing hyw_varid data');
netcdf.putVar(nc,hyw_varid,Mobj.w);

disp('writing zeta_varid data');
netcdf.putVar(nc,zeta_varid,Mobj.surfaceElevation);

% for type 3 nesting need to include weighting
if (Mobj.NEST_TYPE=='TYPE3')

    disp('writing node_weight_varid data');
    netcdf.putVar(nc,node_weight_varid,Mobj.weight_node_val);

    disp('writing element weighting attributes');
    netcdf.putVar(nc,elem_weight_varid,Mobj.weight_elem_val);



end


% close file
netcdf.close(nc);

if(report); fprintf(['end   : ' subname '\n']); end;


end
