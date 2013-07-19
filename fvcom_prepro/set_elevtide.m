function set_elevtide(Mobj,JulianTime,SurfaceElevation,ElevationFile,MyTitle)

% Setup surface elevation tides on the open boundary and dump an
% elevation time series file.
%
% function set_elevtide(Mobj,JulianTime,SurfaceElevation,ElevationFile,MyTitle)
%
% DESCRIPTION:
%    Setup surface elevation tides on the open boundary and dump a NetCDF
%    file.
%
% INPUT
%    Mobj          = Matlab mesh object
%    JulianTime    = Array of time steps in Julian Time (NOT Modified
%    Julian Time)
%    ElevationFile = Output file name
%    MyTitle       = Title in resulting NetCDF file.
%
% OUTPUT:
%
% EXAMPLE USAGE
%    set_elevtide(Mobj,JulianTime,SurfaceElevation,ElevationFile,MyTitle)
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2012-08-08 First version.
%
%==============================================================================
subname = 'set_spectide';
global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

cnt = 0;
ObcNodes = nan(1,sum(Mobj.nObcNodes));
for ob=1:Mobj.nObs
	nObcs = Mobj.nObcNodes(ob);
	for j=1:nObcs
		cnt = cnt + 1;
		ObcNodes(cnt) = Mobj.obc_nodes(ob,j);  % set the open boundary nodes
    end
end

%------------------------------------------------------------------------------
% Dump a surface elevation tide file in NetCDF
%------------------------------------------------------------------------------
write_FVCOM_elevtide(ObcNodes,JulianTime,SurfaceElevation,ElevationFile,MyTitle)

if(ftbverbose); fprintf(['end   : ' subname '\n']); end
