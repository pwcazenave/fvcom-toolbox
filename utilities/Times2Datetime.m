function [ outputDTs ] = Times2Matlab( Times )
%TIMES2DATETIME Converts the 'Times' variable from a FVCOM output file to
%a Matlab Datetime array. This isn't gonna work if the FVCOM model isn't using 
% Julian dates.
%   Input: 'Times' variable, already read from a FVCOM netCDF output file.
%   Output: Array of Matlab Datetime objects, one for each timestep.

% NB time zones are not considered.

%   Simon Waldman / Marine Scotland Science 2018.

global ftbverbose;
if ftbverbose
    [~, subname] = fileparts(mfilename('fullpath'));
    fprintf('\nbegin : %s\n', subname)
end

assert(ischar(Times), 'Input should be a character array.');

% FVCOM can use two different formats for the Times variable, depending on
% which file it's in. So we have to mess around with accounting for both.
% There's probably a better way to do this that doesn't rely on catching
% errors...

try
    outputDTs = datetime( Times', 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS'); 
    if ftbverbose
        disp('Successfully converted Times to datetime using ''yyyy-MM-dd''T''HH:mm:ss.SSSSSS'' format.');
    end
    return;
catch
    if ftbverbose
        disp('Unable to convert the text to datetime using the format ''yyyy-MM-dd''T''HH:mm:ss.SSSSSS''. Trying other format.');
    end
end

try
    outputDTs = datetime( Times', 'InputFormat', 'yyyy/MM/dd'' ''HH:mm:ss.SSSSSS');
    if ftbverbose
        disp('Successfully converted Times to datetime using ''yyyy/MM/dd'' ''HH:mm:ss.SSSSSS'' format.');
    end;
    return;
catch
    if ftbverbose
        disp('Unable to convert the text to datetime using the format ''yyyy/MM/dd'' ''HH:mm:ss.SSSSSS''.');
    end
    error('Failed to convert Times to datetime using either FVCOM format.');
end

error('We should never get to here. Something went wrong.');

end

