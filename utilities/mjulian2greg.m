function [year,month,day,hour,minu,sec,dayweek,dategreg] = mjulian2greg(MJD)
% Convert Modified Julian dates to Gregorian dates.
%
% [year,month,day,hour,minu,sec,dayweek,dategreg] = julian2greg(MJD)
%
% DESCRIPTION:
%   Convert Modified Julian date to Gregorian dates.
%
% INPUT:
%   MJD = list of Modified Julian Dates
%
% OUTPUT:
%   year = Gregorian year
%   month = Gregorian month
%   day = Gregorian day
%   hour = Gregorian hour
%   minu = Gregorian minute
%   sec = Gregorian second
%   dayweek = day of the week for the returned date
%   dategreg = vector of Gregorian dates in the format (D, M, Y, H, M S)
%
% EXAMPLE USAGE
% [year,month,day,hour,minu,sec,dayweek,dategreg] = mjulian2greg(MJD)
%
% Author(s):
%   Geoff Cowles (University of Massachusetts Dartmouth)
%   Karen Amoudry (National Oceanography Centre, Liverpool)
%
% Revision history:
%   2013-09-03 (KJA) Major updates to original version (kept intact in
%   comments at the bottom of this version). The original version
%   frequently generated incorrect times due to rounding errors, with time
%   vectors containing 'impossible' times such as [2008,2,1,0,59,60], or
%   incorrectly rounded times such as [2008,2,1,2,0,1]. The 'impossible'
%   times occasionally caused crashes when used in FVCOM input files. This
%   updated version uses Matlab's inhouse date manipulation functions by
%   accounting for the offset between the origins of Modified Julian Time
%   (1858 11 17 00:00:00) and Matlab datenum time (0000 01 01 00:00:00).
%
%==========================================================================

global ftbverbose
report = false;
if(ftbverbose); report = true; end
subname = 'mjulian2greg';
if(report); fprintf('\n'); end
if(report); fprintf(['begin : ' subname '\n']); end

% Calculate the offset between the origins of Modified Julian Time and
% Matlab datenum time.
MJD_base=datenum(1858,11,17,0,0,0);

% Convert from MJD to Gregorian time
[year,month,day,hour,minu,sec] = datevec(MJD+MJD_base);

% Calculate the day of the week
[~,dayweek] = weekday(MJD+MJD_base,'long');

% Convert to dategreg format (D, M, Y, H, M S)
dategreg = [day, month, year, hour, minu, sec];

if(report); fprintf(['end   : ' subname '\n']); end;

%% Original version of mjulian2greg.m
%This function converts Modified Julian dates to Gregorian dates.

%[year,month,day,hour,minu,sec,dayweek,dategreg] = julian2greg(MJD+2400000.5);
