function swan2netcdf(matfile,ncfile,basename,first_time,last_time,increment,isforcing,kn)
%matfile = 'c4.mat';
%ncfile = 'case_s0004.nc';
%basename = 'skg4.3';
%first_time = '20090618_000000';
%last_time = '20090621_000000';
%increment = 3600;
%isforcing = true;
%kn = .003;
% Convert a SWAN output file (Matlab) into a NetCDF file
%
% function swan2netcdf(matfile,ncfile,basename,first_time,last_time,increment);
%
% DESCRIPTION:
%    read output from unstructured SWAN model (currently 40.82) and
%    dump to a NetCDF file which is far more useful than a Matlab file.
%
% INPUT
%   matfile  = Unstructured SWAN Matlab file
%   ncfile   = NetCDF file for output
%   basename = prefix for SWAN mesh, bathymetry, connectivity files
%   first_time:  first time frame in Matlab object frame time
%   last_time:   last time frame in Matlab object frame time
%   increment:   increment in seconds
%   isforcing:   converts NaNs from .mat file to 1.0's or 0.0's
%   kn:          Nikuradse roughness in meters  [OPTIONAL, default = .01]
% OUTPUT:
%    NetCDF file containing:
%      a.) time in modified Julian day
%      b.) significant wave height (hs)
%      c.) wave direction
%      d.) mesh
%      e.) peak period
%      f.) U10
%      g.) V10
%      h.) bottom orbital velocity
%      i.) bottom period
%      j.) H/h ratio
%      k.) wave-induced bed stress
%
% EXAMPLE USAGE
%   swan2netcdf('gom1.mat','gom1.nc','gom1','20070101_000000','20070131_000000',3600,true,.045)
%     this converts gom1.mat to gom1.nc using SWAN grid files gom1.ele, gom1.bot
%     and gom1.node from Jan 1, 2007 00:00:00 to Jan 31, 2007 00:00:00 in increments
%     of 1 hour.
%
% NOTE
%    routine is not refined, e.g. will not check if files exist and will
%    probably crash if you do not have the variables above in the SWAN
%    output file.
%    You will need approximately the following BLOCK command in your SWAN runfile
%
%    BLOCK 'COMPGRID' NOHEAD 'gom1.mat' LAY 3 XP YP DEP HS RTP TPS DIR WLEN &
%                WINDV OUTPUT 20070101_000000 3600 SEC
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Eric Holmes (University of Massachusetts Dartmouth)
%
% Revision history
%    09-29-2010 -Now warns you if variable was not found in SWAN output
%    .mat file and added two variables from SWAN output Ubot and TmBot.
%    12-01-2010 Added 'isforcing' option
%==============================================================================

% set Nikaradse roughness to 1 cm if user does not specify
if(~exist('kn'))
  kn = .01;
end;
z0 = kn/30.;  %hydraulic roughness

eval(['load ' matfile]);              % load binary file containing SWAN results
% obtained using BLOCK command with COMPGRID-set
% load the connectivity
elefile=[basename '.ele'];
fid = fopen(elefile);                 % load TRIANGLE element based connectivity file
[nele] = fscanf(fid,'%i',[1 3]);     % get number of triangles
jnk = fscanf(fid,'%i',[4 nele(1)])'; % get connectivity table
tri = jnk(:,2:4);

% load the bathymetry
bathfile=[basename '.bot'];
h = textread(bathfile,'%f\n');
node = prod(size(h));

% check if variable range exists
f = first_time;
l = last_time;
tbeg = datenum(str2num(f(1:4)),str2num(f(5:6)),str2num(f(7:8)),str2num(f(10:11)),str2num(f(12:13)),str2num(f(14:15)));
tend = datenum(str2num(l(1:4)),str2num(l(5:6)),str2num(l(7:8)),str2num(l(10:11)),str2num(l(12:13)),str2num(l(14:15)));
tinc  = increment/(24*3600);
times = tbeg:tinc:tend;
ntimes = prod(size(times));
for i=1:ntimes;
    %  datestr(times(i))
    date = datestr(times(i),30);
    vname = ['Hsig_',date(1:8),'_',date(10:15)];
    if(exist(vname) == 0)
        error('variable frame %s\n does not exist',vname)
    end;
end;


% dump the netcdf file
nc = netcdf(ncfile, 'clobber');
nc.info = 'created from SWAN unstructured grid output file';
nc.source = 'fvcom grid (unstructured) surface forcing';

% dimensions
nc('nele') = nele;
nc('node') = node;
nc('three') = 3;
nc('time') = 0;

% variables
date = datestr(times(1),30);



nc{'x'} = ncfloat('node');
nc{'x'}.long_name = 'nodal x-coordinate';
nc{'x'}.units     = 'm';

nc{'y'} = ncfloat('node');
nc{'y'}.long_name = 'nodal y-coordinate';
nc{'y'}.units     = 'm';

nc{'h'} = ncfloat('node');
nc{'h'}.long_name = 'Bathymetry';
nc{'h'}.units     = 'm';

nc{'d'} = ncfloat('time','node');
nc{'d'}.long_name = 'Depth';
nc{'d'}.units     = 'm';

nc{'nv'} = ncint('three','nele');
nc{'nv'}.long_name = 'nodes surrounding element';


nc{'time'} = ncfloat('time');
nc{'time'}.long_name = 'time';
nc{'time'}.units = 'days since 1858-11-17 00:00:00';
nc{'time'}.format = 'modified julian day (MJD)';
nc{'time'}.time_zone = 'UTC';

nc{'Itime'} = ncint('time');
nc{'Itime'}.units = 'days since 1858-11-17 00:00:00';
nc{'Itime'}.format = 'modified julian day (MJD)';
nc{'Itime'}.time_zone = 'UTC';

nc{'Itime2'} = ncint('time');
nc{'Itime2'}.units = 'msec since 00:00:00';
nc{'Itime2'}.time_zone = 'UTC';

vname = ['Hsig_',date(1:8),'_',date(10:15)];
if(exist(vname) == 0)
    warning('Variable Hsig does not exist not being added to netCDF output')
else
    nc{'hs'} = ncfloat('time','node');
    nc{'hs'}.long_name = 'Significant Wave Height';
    nc{'hs'}.units     = 'm';
end;

vname = ['Dir_',date(1:8),'_',date(10:15)];
if(exist(vname) == 0)
    warning('Variable Dir does not exist not being added to netCDF output')
else
    nc{'wdir'} = ncfloat('time','node');
    nc{'wdir'}.long_name = 'Wave  Direction';
    nc{'wdir'}.units     = 'degree';
end;

vname = ['RTpeak_',date(1:8),'_',date(10:15)];
if(exist(vname) == 0)
    warning('Variable RTpeak does not exist not being added to netCDF output')
else
    nc{'tpeak'} = ncfloat('time','node');
    nc{'tpeak'}.long_name = 'Relative Peak Period';
    nc{'tpeak'}.units     = 's';
end;

vname = ['Windv_x_',date(1:8),'_',date(10:15)];
if(exist(vname) == 0)
    warning('Variable WindV_x_ does not exist not being added to netCDF output')
else
    nc{'U10'} = ncfloat('time','node');
    nc{'U10'}.long_name = 'Wind Velocity x-direction';
    nc{'U10'}.units     = 'm/s';
end;

vname = ['Windv_y_',date(1:8),'_',date(10:15)];
if(exist(vname) == 0)
    warning('Variable WindV_y_ does not exist not being added to netCDF output')
else
    nc{'V10'} = ncfloat('time','node');
    nc{'V10'}.long_name = 'Wind Velocity y-direction';
    nc{'V10'}.units     = 'm/s';
end;

vname = ['Wlen_',date(1:8),'_',date(10:15)];
if(exist(vname) == 0)
    warning('Variable Wlen does not exist not being added to netCDF output')
else
    nc{'wlen'} = ncfloat('time','node');
    nc{'wlen'}.long_name = 'wavelength';
    nc{'wlen'}.units     = 'm';
end;

vname = ['Ubot_',date(1:8),'_',date(10:15)];
if(exist(vname) == 0)
    warning('Variable Ubot does not exist not being added to netCDF output')
else
    nc{'Ubot'} = ncfloat('time','node');
    nc{'Ubot'}.long_name = 'Bottom Orbital Velocity';
    nc{'Ubot'}.units     = 'm/s';
end;


vname = ['TmBot_',date(1:8),'_',date(10:15)];
if(exist(vname) == 0)
    warning('Variable TmBot does not exist not being added to netCDF output')
else
    nc{'TmBot'} = ncfloat('time','node');
    nc{'TmBot'}.long_name = 'Bottom Wave Period';
    nc{'TmBot'}.units     = 's';
end;

vname1 = ['Hsig_',date(1:8),'_',date(10:15)];
vname2 = ['Depth_',date(1:8),'_',date(10:15)];
if(exist(vname1) == 0 | exist(vname2)==0)
    warning('Bot Hsig and Depth do not exist so will not dump H/h')   
else
    nc{'H_over_h'} = ncfloat('time','node');
    nc{'H_over_h'}.long_name = 'Ratio of Hsig to depth';
    nc{'H_over_h'}.units     = '-';
end;

vname1 = ['TmBot_',date(1:8),'_',date(10:15)];
vname2 = ['Ubot_',date(1:8),'_',date(10:15)];
if(exist(vname1) == 0 | exist(vname2)==0)
    warning('Bot TmBot and Ubot do not exist so will not dump bed stress')
else
    nc{'tau_w'} = ncfloat('time','node');
    nc{'tau_w'}.long_name = 'wave-induced bed stress'; 
    nc{'tau_w'}.units     = 'm^2/s^2';
end;

% static vars
nc{'x'}(:) = Xp;
nc{'y'}(:) = Yp;
nc{'h'}(:) = h;
nc{'nv'}(:,:) = tri';

badpts = [];
% dump dynamic vars
for i=1:ntimes;
    
    fprintf('processing time %s\n',datestr(times(i)));
    
    if(isforcing==1)
	    date = datestr(times(i),30);
		vname = ['Hsig_',date(1:8),'_',date(10:15)]; var1 = eval(vname)';
		vname = ['RTpeak_',date(1:8),'_',date(10:15)];var2 = eval(vname)';
		var3  = var1+var2;
		badpts = find(isnan(var3));
		clear var1;
		clear var2;
		clear var3;
	end; 
	
    %time
    shift = 678942.;  % datenum(2010,1,1,0,0,0)-greg2mjulian(2010,1,1,0,0,0);
    time  = times(i) - shift;
    nc{'time'}(i) = time;
    nc{'Itime'}(i) = floor(time);
    nc{'Itime2'}(i) = mod(time,1)*24*3600*1000.;
    
    
    %hs
    date = datestr(times(i),30);
    vname = ['Hsig_',date(1:8),'_',date(10:15)];
    if(exist(vname) == 0)
        %fprintf('variable frame %s\n does not exist',vname)
    else
        var = eval(vname)';
        if (isforcing == 1); var(badpts) = 0.; end;
        nc{'hs'}(i,:) = var;
    end;
    
    %tp
    date = datestr(times(i),30);
    vname = ['RTpeak_',date(1:8),'_',date(10:15)];
    if(exist(vname) == 0)
        %fprintf('variable frame %s\n does not exist',vname)
    else
        var = eval(vname)';
		if (isforcing == 1);  var(badpts) = 1.; end;
        nc{'tpeak'}(i,:) = var;
    end;
    
    
    %depth
    date = datestr(times(i),30);
    vname = ['Depth_',date(1:8),'_',date(10:15)];
    if(exist(vname) == 0)
        %fprintf('variable frame %s\n does not exist',vname)
    else
        nc{'d'}(i,:) = eval(vname)';
    end;
    
    % wave dir
    date = datestr(times(i),30);
    vname = ['Dir_',date(1:8),'_',date(10:15)];
    if(exist(vname) == 0)
        %fprintf('variable frame %s\n does not exist',vname)
    else
        var = eval(vname)';
        if (isforcing == 1); var(badpts) = 0.; end;
        nc{'wdir'}(i,:) = var;
    end;
    
    % U10
    date = datestr(times(i),30);
    vname = ['Windv_x_',date(1:8),'_',date(10:15)];
    if(exist(vname) == 0)
        %fprintf('variable frame %s\n does not exist',vname)
    else
        var = eval(vname)';
        nc{'U10'}(i,:) = var;
    end;
    
    % V10
    date = datestr(times(i),30);
    vname = ['Windv_y_',date(1:8),'_',date(10:15)];
    if(exist(vname) == 0)
        %fprintf('variable frame %s\n does not exist',vname)
    else
        var = eval(vname)';
        if (isforcing == 1)
            var(isnan(var)) = 0.0;
        end
        nc{'V10'}(i,:) = var';
    end;
    
    % orbital velocity
    date = datestr(times(i),30);
    vname = ['Ubot_',date(1:8),'_',date(10:15)];
    if(exist(vname) == 0)
        %fprintf('variable frame %s\n does not exist',vname)
    else
        var = eval(vname)';
        if (isforcing == 1); var(badpts) = 0.; end;
        nc{'Ubot'}(i,:) = var';
    end;
    
    % wavelength
    date = datestr(times(i),30);
    vname = ['Wlen_',date(1:8),'_',date(10:15)];
    if(exist(vname) == 0)
        %fprintf('variable frame %s\n does not exist',vname)
    else
        var = eval(vname)';
        if (isforcing == 1); var(badpts) = 1.; end;
        nc{'wlen'}(i,:) = var';
    end;
    
    % bottom wave period
    date = datestr(times(i),30);
    vname = ['TmBot_',date(1:8),'_',date(10:15)];
    if(exist(vname) == 0)
        %fprintf('variable frame %s\n does not exist',vname)
    else
        var = eval(vname)';
        if (isforcing == 1); var(badpts) = 1.; end;
        nc{'TmBot'}(i,:) = var';
    end;

    % hsig/depth
    date = datestr(times(i),30);
    vname1 = ['Hsig_',date(1:8),'_',date(10:15)];
    vname2 = ['Depth_',date(1:8),'_',date(10:15)];
    if(exist(vname1) == 0 | exist(vname2)==0)
        %fprintf('variable frame %s\n does not exist',vname)
    else
        var1 = eval(vname1)';
        var2 = eval(vname2)';
        var  = var1./var2;
        if (isforcing == 1); var(badpts) = 0.; end;
        nc{'H_over_h'}(i,:) = var';
    end;
    
    % wave-induced bed stress  
    date = datestr(times(i),30);
    vname1 = ['TmBot_',date(1:8),'_',date(10:15)];
    vname2 = ['Ubot_',date(1:8),'_',date(10:15)];
    if(exist(vname1) == 0 | exist(vname2)==0)
        %fprintf('variable frame %s\n does not exist',vname)
    else
        %compute wave induced bed stress using Soulsby formulas
        %friction factor (fw) is eq 62A 
        %tau_w = 0.5*fw*Ubot^2
        var1 = eval(vname1)';
        var2 = eval(vname2)';
        omega_wave = (2.0*pi./max(var1,.05))'; %wave-orbital frequency
        var=0.5*1.39*((omega_wave*z0).^.52)'.*(var2.^(2.0-.52));
        if (isforcing == 1); var(badpts) = 0.; end;
        nc{'tau_w'}(i,:) = var';
    end;
    %error('hog')

end;


nc = close(nc);
