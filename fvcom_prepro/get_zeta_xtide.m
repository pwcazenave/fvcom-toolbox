function [times,zeta] = get_zeta_xtide(Xtide_Location,tbeg,tend)

% Extract an elevation time series from a location using Xtide 
%
% function [time,zeta] = gen_zeta_xtide(Xtide_Location,tbeg,tend)  
%
% DESCRIPTION:
%    Extract time and an elevation series from a location using xtide
%
% INPUT:
%   Xtide_Location = Location in Xtide Database 
%   tbeg           = begin time in modified Julian days UTC
%   tend           = end time in modified Julian days UTC
%
% OUTPUT:
%    time  = time sequence for elevation series in modified Julian days, UTC
%    zeta  = elevation series in meters
%
% EXAMPLE USAGE
%   [time,zeta] = get_zeta_xtide('Sandy Point, Whidbey Island',54191,54466) 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================


% conversions for dates and lengths 
unixstart      = greg2mjulian(1970,1,1,0,0,0);  %xtide is based on days from this date
meters2feet    = 3.2808399;

% paths to executable
setpath        = ' export HFILE_PATH=~gcowles/Packages/xtide/harmonics_files; ';
ttide          = ' /usr/local/bin/tide -b "2005-01-01 00:00" -e "2012-01-01 00:00" -l ';
opts           = ' -s "00:06"  -m r -z > ';

%----------------------------------------------------------------------------------
% extract time series from xtide at station 
%----------------------------------------------------------------------------------
cmd = [setpath ttide  '"' char(Xtide_Location) '"' opts  '"tidefile"'];
system(cmd);
[times,zeta] = textread('tidefile','%f %f\n');
system('\rm tidefile');

%-------------------------------------------------------------
% process xtide data  
%    - convert to meters
%    - shift to MSL
%    - shift to Julian Day, UTC
%-------------------------------------------------------------

% convert feet => meters 
zeta = zeta/meters2feet;

% shift the vertical datum to center around mean 
zeta = zeta - mean(zeta);

% convert the time from unix time to modified Julian day UTC
times = unixstart + times/(3600.*24.);

% determine begin/end frames
[mbeg,ibeg] = min(abs(times-tbeg));
[mbeg,iend] = min(abs(times-tend));

times = times(ibeg:iend);
zeta  = zeta(ibeg:iend);



