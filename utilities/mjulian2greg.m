function [year,month,day,hour,minu,sec,dayweek,dategreg] = julian2greg(MJD)
% This function converts Modified Julian dates to Gregorian dates.

[year,month,day,hour,minu,sec,dayweek,dategreg] = julian2greg(MJD+2400000.5);
