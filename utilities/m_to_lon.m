function dlon = m_to_lon(dx, alat)
%
% dlon = longitude difference in degrees
% dx   = longitude difference in meters
% alat = average latitude between the two fixes

rlat = alat * pi/180;
p = 111415.13 * cos(rlat) - 94.55 * cos(3 * rlat);
dlon = dx ./ p;
