function all_data = concat_by_struct(in_struct,last)

% A function to concatenate data contained within a structure.
%     does something like cat(1,data(:).*);  
% function all_data = concat_by_struct(in_struct,last)
% 
% DESCRIPTION:
%    Concatenates data contained within a structure. It only works on the first dimension unless
%     logical variable "last" is set to true.
%
% INPUT 
%   in_struct  = Structured variable with dimension such as 
%		in_struct(1:10).data1, in_struct(1:10).data2...
%   last       = logical variable. Enable to work on last dimension only
%
% OUTPUT:
%    all_data = struct with first dimension removed.
%
% EXAMPLE USAGE
%    all_data = concat_by_struct(in_struct,true)
%
% Author(s):
%    Ricardo Torres (Plymouth Marine Laboratory) 
%
% Revision history
%
%   2018-03-22 First version in fvcom-toolbox.
%
%==============================================================================

subname = mfilename;
global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

vars = fieldnames(in_struct);
maxl=1;sizes=1;
for kk = 1:size(vars,1)
% find largest record
eval(['xyz=size(in_struct(1).',char(vars(kk)) ');'])
sizes=[sizes;xyz(:)];;
end

for kk = 1:size(vars,1)
    if eval(['isstruct(in_struct(1).' char(vars(kk)) ')'])
        in_vars = fieldnames(in_struct);
        eval(['all_data.' char(vars(kk)) ' = concat_by_struct([in_struct.' char(vars(kk)) '],last);']);
    else
        eval(['xyz=ndims(in_struct(1).',char(vars(kk)) ');'])
        if last % If concatenating on the last dimension
            
                eval(['all_data.' char(vars(kk)) ' = cat(xyz,in_struct(:).' char(vars(kk)) ');']);
        else % concatenate on the first dimension
                    eval(['all_data.' char(vars(kk)) ' = cat(1,in_struct(:).' char(vars(kk)) ');']);
        end
    end
end
    
if ftbverbose
  fprintf('end   : %s\n', subname)
end
