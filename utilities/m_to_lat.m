function dlat = m_to_lat(dy, alat)
%
% dy   = latitude difference in meters
% dlat = latidute difference in degrees
% alat = average latitude between the two fixes
% Reference: American Practical Navigator, Vol II, 1975 Edition, p 5 

rlat = alat * pi/180;
m = 111132.09 * ones(size(rlat)) - ...
    566.05 * cos(2 * rlat) + 1.2 * cos(4 * rlat);
dlat = dy ./ m;
