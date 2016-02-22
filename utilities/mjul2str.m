function strout = mjul2str(MJD,noyear)
% Convert a modified Julian day to a Matlab datestr style string 

MJD = double(MJD);

mjul2matlab = 678942; %difference between modified Julian day 0 and Matlab day 0
if(~exist('noyear'))
strout = datestr(MJD+mjul2matlab);
else
strout = datestr(MJD+mjul2matlab);
strout = strout(1:end-5)
end;
