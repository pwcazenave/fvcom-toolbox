function scoord=read_scoord_params(filename)
% function to read hc, cc, theta and bb from scoord_params.dat file
% filename='/users/modellers/rito/Models/MEDINA/polcoms/scoord_params.dat'
fid = fopen(filename);
C = textscan(fid, '%f%*[^\n]');
fclose(fid);
scoord.hc=C{1}(1);
scoord.cc=C{1}(2);
scoord.theta=C{1}(3);
scoord.bb=C{1}(4);
return
