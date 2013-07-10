function [mjulianday] =greg2mjulian(year,month,day,hour,mint,sec)
% mjulianday = greg2mjulian(yyyy,mm,dd,HH,MM,SS)
%
% Convert the Gregorian dates to Modified Julian dates.
%
% See GREG2JULIAN.

mjulianday = greg2julian(year,month,day,hour,mint,sec)-2400000.5;
