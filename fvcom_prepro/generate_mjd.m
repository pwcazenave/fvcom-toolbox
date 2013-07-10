%--------------------------------------------------------------
% dump to netcdf file
%--------------------------------------------------------------

start_day = greg2mjulian(2007,4,1,1,0,0);
end_day   = greg2mjulian(2007,7,1,0,0,0);

time = start_day:(1./24):end_day;

% open boundary forcing
nc = netcdf('gom1v10_decelles_2007_time.nc', 'clobber');       


% dimensions
nc('time') = 0; 

% variables
nc{'time'} = ncfloat('time');
nc{'time'}.long_name = 'time';  
nc{'time'}.units     = 'days since 0.0';  
nc{'time'}.time_zone = 'none';  

nc{'Itime'} = ncint('time');
nc{'Itime'}.units     = 'days since 0.0';  
nc{'Itime'}.time_zone = 'none';  

nc{'Itime2'} = ncint('time');
nc{'Itime2'}.units     = 'msec since 00:00:00';
nc{'Itime2'}.time_zone = 'none';  


% dump dynamic data
ntimes = numel(time);
nc{'time'}(1:ntimes) = time(1:ntimes);
nc{'Itime'}(1:ntimes) = floor(time(1:ntimes));
nc{'Itime2'}(1:ntimes) = mod(time(1:ntimes),1)*24*3600*1000.;

nc = close(nc);    


