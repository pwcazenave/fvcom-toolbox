function [varlist] = read_fabm_variables(conf)
% For a given configuration, read the ersem variables required to run the
% model
%
% [varlist] = get_fabm_variables(conf)
%
% DESCRIPTION:
%   The conf struct contains the parameters needed to locate the restart
%   file. We only need the variable names...
%
% INPUT:
%   conf - struct. With fields:
%       * restart_file - FVCOM restart file with FABM variables (i.e. a donor file).
% OUTPUT:
%   varlist - cell array of ERSEM variables.
%
% EXAMPLE USAGE:
%   [varlist] = get_fabm_variables(conf);
%
% Author(s):
%   Ricardo torres (Plymouth Marine Laboratory)
%
% Revision history:
%   2017-03-10 - First version 
%
%==========================================================================
% hard coding these is a bad habit!! These should be identical to the ones
% in get_POLCOMS_ERSEM.m
% conf.casename = 'aqua_v14';
% conf.base_dir = '/users/modellers/rito';
% conf.project_dir = fullfile(conf.base_dir, '/Models/git/fvcom/rosa-implementation/');
% conf.restart_file = fullfile(conf.project_dir, 'matlab/restarts/', ...
%     sprintf('%s_ersem_restart_donor.nc', conf.casename));


subname = 'read_fabm_variables';
global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

% check file exist
% 

if exist(conf.restart_file,'file')
% read ersem variables in donor restart file
info_donor=ncinfo(conf.restart_file);
else
    warning(['File ',conf.restart_file,' Not found. Returning'])
    return
end
% Generally the last fvcom variable although this could change with a 
% different configuration of output variables. Better check on a case by case 
% There is no FABM reference in the variable attributes to possibly do it differently
last_fvcom = 'wet_cells_prev_ext';
idx = find(strcmpi({info_donor.Variables.Name},last_fvcom));
ii=1;varlist={};
for ff=idx+1:length(info_donor.Variables)
    % brute force search...
    res = strfind(info_donor.Variables(ff).Name,['_']);
    if ~isempty(res)
        % variable seems to have the construct of FABM ...
        % check is not a benthic variable
        if (length(info_donor.Variables(ff).Dimensions)==3)
            varlist(ii)={info_donor.Variables(ff).Name};
            ii = ii+1;
            
        end
    end
end

return
