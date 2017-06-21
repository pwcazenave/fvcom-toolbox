function [h] = read_fvcom_bath(bathfile) 

% Read fvcom bathymetry file 
%
% [h] = function read_fvcom_bath(bathfile)
%
% DESCRIPTION:
%    Read FVCOM Bathymetry file 
%
% INPUT [keyword pairs]:  
%   'bathfile'  = fvcom bathymetry file
%
% OUTPUT:
%    h = bathymetry vector
%
% EXAMPLE USAGE
%    Mobj = read_fvcom_bath('tst_dep.dat')
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Rory O'Hara Murray (Marine Scotland Science)
%
% Revision history
%    2014-11-19 Remove loops to speed up reading in the file.
%   
%==============================================================================

global ftbverbose
[~, subname] = fileparts(mfilename('fullpath'));
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

%------------------------------------------------------------------------------
% read in the FVCOM bathymetry data
%------------------------------------------------------------------------------
fid = fopen(bathfile,'r');
if(fid  < 0)
  error(['file: ' bathfile ' does not exist']);
end;
C = textscan(fid, '%s %s %s %d', 1);
Nverts = C{4};
h = zeros(Nverts,1);
fprintf('reading bathymetry file\n');
fprintf('# nodes %d\n',Nverts);
C = textscan(fid,' %f %f %f',Nverts);
h = C{3};
fprintf('min depth %f max depth %f\n',min(h),max(h));
fprintf('bathymetry reading complete\n');
fclose(fid);



if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end;


