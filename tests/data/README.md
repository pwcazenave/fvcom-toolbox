Notes
=====

To extract only the required variables from the full POLCOMS file, use the NetCDF Kitchen Sink toolkit as follows:

% ncks -a -v time,pdepthD,ETWD,x1XD,lon,lat <infile> <outfile>
