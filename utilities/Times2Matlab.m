function [ sdTimes ] = Times2Matlab( Times )
%TIMES2MATLAB Converts the 'Times' variable from a FVCOM output file to
%Matlab serial dates. This isn't gonna work if the FVCOM model isn't using 
% Julian dates.
%   Input: 'Times' variable, already read from a FVCOM netCDF output file.
%   Output: Array of Matlab serial date numbers, one for each timestep.

%   Simon Waldman / Heriot-Watt University 2016.

sdTimes = datenum(Times(1:23,:)', 'yyyy-mm-ddTHH:MM:SS.FFF'); 

% Note - FVCOM includes microseconds, but datenum doesn't recognise that
% level of precision. To be correct I should round them, but instead I
% truncate them because it's simple. If you care about +/- 1 millisecond in your model results,
% don't use this function.

end

