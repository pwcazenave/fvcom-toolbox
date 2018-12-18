function [ outputDTs ] = Times2Matlab( Times )
%TIMES2DATETIME Converts the 'Times' variable from a FVCOM output file to
%a Matlab Datetime array. This isn't gonna work if the FVCOM model isn't using 
% Julian dates.
%   Input: 'Times' variable, already read from a FVCOM netCDF output file.
%   Output: Array of Matlab Datetime objects, one for each timestep.

% NB time zones are not considered.

%   Simon Waldman / Marine Scotland Science 2018.

assert(ischar(Times), 'Input should be a character array.');

outputDTs = datetime( Times', 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS'); 

end

