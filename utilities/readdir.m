function [varargout] = readdir();
%Read working directory with no arguments. Output depends on the function call
% 
%
fdir = fopen('./dirnames.dir','r');
%strippintg comment lines if any
header=fgetl(fdir);
if header(1)=='%'
   while header(1)=='%'
      header=fgetl(fdir);
   end
   varargout(1)= {header};
   for i=2:nargout 
      varargout(i)= {fgetl(fdir)};
   end
else
   frewind(fdir)
   for i=1:nargout 
      varargout(i)= {fgetl(fdir)};
   end
end
fclose(fdir);
return