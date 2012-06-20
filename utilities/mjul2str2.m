function strout = mjul2str2(MJD)
% Convert a modified Julian day to a Matlab datestr style string 
% Do not return year

mjul2matlab = 678942; %difference between modified Julian day 0 and Matlab day 0
strout = datestr(MJD+mjul2matlab);
strout = strout(1:end-5);
