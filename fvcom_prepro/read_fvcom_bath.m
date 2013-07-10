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
%
% Revision history
%   
%==============================================================================

subname = 'read_fvcom_bath';
global ftbverbose
if(ftbverbose)
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

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
for i=1:Nverts
  C = textscan(fid, '%f %f %f', 1);
  h(i) = C{3};
end;
fprintf('min depth %f max depth %f\n',min(h),max(h));
fprintf('bathymetry reading complete\n');
fclose(fid);



if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end;


