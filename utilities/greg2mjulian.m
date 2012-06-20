function [mjulianday] =greg2mjulian(year,month,day,hour,mint,sec)
% This function converts the Gregorian dates to Modified Julian dates.
% -----------------------------------------------------------------------------------------------------------

mjulianday = greg2julian(year,month,day,hour,mint,sec)-2400000.5;
