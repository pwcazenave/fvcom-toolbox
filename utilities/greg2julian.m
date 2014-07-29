function [JD,julianday] =greg2julian(year,month,day,hour,min,sec)
% This function converts the Gregorian dates to Julian dates.
%
% 0. Syntax:
% [JD,julianday] = juliandate(year,month,day,hour,min,sec)
%
% 1. Inputs:
%     year, month, day = date in Gregorian calendar.
%     hour,min,sec = time at universal time.
%
% 2. Outputs:
%     JD = Julian date.
%     julianday = day of week.
%
% 3. Example:
%  >> [a,b] = greg2julian(2006,5,30,2,30,28)
%  a =
%
%           2453885.60449074
%  b =
%
%  Tuesday
%
%  4. Notes:
%     - For all common era (CE) dates in the Gregorian calendar, and for more
%     information, check the referents.
%     - The function was tested, using  the julian date converter of U.S. Naval Observatory and
%     the results were similar. You can check it.
%     - Trying to do the life... more easy with the conversions.
%
% 5. Referents:
%     Astronomical Applications Department. "Julian Date Converter". From U.S. Naval Observatory.
%               http://aa.usno.navy.mil/data/docs/JulianDate.html
%     Duffett-Smith, P. (1992).  Practical Astronomy with Your Calculator.
%               Cambridge University Press, England:  pp. 9.
%     Seidelmann, P. K. (1992). Explanatory Supplement to the Astronomical Almanac.
%               University Science Books, USA.  pp. 55-56.
%      Weisstein, Eric W.  "Julian Date".  From World of Astronomy--A Wolfram Web Resource.
%               http://scienceworld.wolfram.com/astronomy/JulianDate.html
%
% Gabriel Ruiz Mtz.
% May-2006
%
% Modifications:
% 04/06/06: To find the days, it was only changed the loop to a cell array. Thanks to Jerome.
% 2014-07-29 Fix input arguments check (Pierre Cazenave, Plymouth Marine
%   Laboratory).
% ------------------------------------------------------------------------------------------------------------

   narginchk(6,6)
   timeut = hour + ( min / 60 ) + ( sec / 3600 );

   %For common era (CE), anno domini (AD)
   JD = ( 367 * year ) - floor ( 7 * ( year + floor( ( month + 9 ) / 12 ) ) / 4 ) - ...
                     floor( 3 * ( floor( ( year + ( month - 9 ) / 7 ) / 100 ) + 1 ) / 4 ) + ...
                     floor( ( 275 * month ) / 9 ) + day + 1721028.5 + ( timeut / 24 );
   a = ( JD + 1.5 ) / 7;
   frac = a - floor(a);
   n = floor(frac * 7) ;
   julianday ={ 'Sunday' 'Monday' 'Tuesday' 'Wednesday' 'Thursday' 'Friday' 'Saturday'};
   julianday = julianday{n+1};


