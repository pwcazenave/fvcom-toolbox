function [Mobj]  = write_FVCOM_stations(Mobj,filename)

% Add a set of stations at which FVCOM will output time series. 
%
% function write_FVCOM_stations(Mobj,filename)
%
% DESCRIPTION:
%    Given a mesh object with time series positions and names
%    (Mobj.Position and Mobj.Names from add_stations_list.m), write out to
%    ASCII file filename.
%
% INPUT
%    Mobj = Matlab mesh object
%    filename = FVCOM stations file name
%
% OUTPUT:
%    FVCOM stations file: filename
%
% EXAMPLE USAGE
%    write_FVCOM_stations(Mobj, 'tst_stations.dat')
%
% Author(s):  
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2012-11-30 First version.
%   
%==========================================================================
subname = 'write_FVCOM_stations';
global ftbverbose
if(ftbverbose)
  fprintf('\n'); fprintf(['begin : ' subname '\n']);
end

%--------------------------------------------------------------------------
% Parse input arguments
%--------------------------------------------------------------------------
if exist('Mobj', 'var') ~= 1 || exist('filename', 'var') ~= 1
	error('arguments to write_FVCOM_grid are incorrect')
end

%--------------------------------------------------------------------------
% Dump the file
%--------------------------------------------------------------------------
if strcmpi(Mobj.nativeCoords, 'cartesian')
	x = Mobj.Positions(:,3);
	y = Mobj.Positions(:,4);
elseif strcmpi(Mobj.nativeCoords, 'spherical')
	x = Mobj.Positions(:,1);
	y = Mobj.Positions(:,2);
else
    error('Unknown native coordinate system string: %s', Mobj.nativeCoords)
end

if(ftbverbose)
    fprintf('writing FVCOM gridfile %s\n',filename);
end

fid = fopen(filename,'w');
fprintf(fid, ' No           X        Y      Node (Cell)        Station Name\n');

for s=1:length(Mobj.stations)
    fprintf(fid, '%i %f %f %i %f %s\n', cell2mat(Mobj.stations{s}(1:5)), char(Mobj.stations{s}(6)));
end   

fclose(fid);

if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end

