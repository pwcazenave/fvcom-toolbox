function [year,month,day,hour,minu,sec,dayweek,dategreg] = mjulian2greg(MJD)
% This function converts Modified Julian dates to Gregorian dates.
%
% [yr, mon, day, hr, mins, sec, dayweek, dategreg] = mjulian2greg(MJD);

[year,month,day,hour,minu,sec,dayweek,dategreg] = julian2greg(MJD+2400000.5);
