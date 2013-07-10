function dx = lon_to_m(dlon, alat)
% dx = lon_to_m(dlon, alat)
% dx   = longitude difference in meters
% dlon = longitude difference in degrees
% alat = average latitude between the two fixes

rlat = alat * pi/180;
p = 111415.13 * cos(rlat) - 94.55 * cos(3 * rlat);
dx = dlon .* p;
