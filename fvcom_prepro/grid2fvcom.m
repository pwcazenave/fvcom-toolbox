function fvcom = grid2fvcom(Mobj,vars,wind)
% Interpolate regularly gridded wind speed data onto a given FVCOM grid
%
% grid2fvcom(Mobj,vars,wind,fvcom_forcing_file,infos)
% 
% DESCRIPTION:
%   Takes a given NCEP reanalysis grid file and interpolates the U10 and
%   V10 values onto the specified FVCOM grid file. 
%   
% INPUT:
%   Mobj - MATLAB mesh object
%   vars - a cell array of the variables to be interpolated on the FVCOM
%   grid in Mobj (e.g. uwnd, U10, vwnd, V10 etc.).
%   wind - a struct which contains the following arrays:
%       x - x data (probably best in cartesian for the interpolation)
%       y - y data (probably best in cartesian for the interpolation)
%       The struct must also contain all the variables defined in vars.
%       time - time vector (in Modified Julian Days)
% 
% OUTPUT:
%   fvcom - struct of the interpolated data values at the model nodes and
%   element centres.
% 
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
% 
% Revision history:
%   2012-10-15 First version based on ncep2fvcom_U10V10.m in the
%   fvcom-toolbox.
%   2012-10-16 Removed the code to read the NCEP file. Instead, farmed that
%   out to a new function (read_NCEP_wind) so that the relevant section can
%   be more readily extracted (rather than using the entire globe's data:
%   it's easier to subsample and provide the subsampled data here). 
%   2012-10-17 Add outputs to the function for use in visualisation. 
%   2012-10-19 Add wind struct as input rather than separate u, v, time and
%   lat/long arrays. Makes invocation a bit cleaner.
%   2012-11-01 Farmed out the creation of the NetCDF file to
%   write_FVCOM_forcing.m and made this purely an interpolation script. 
% 
%==========================================================================

warning off

if nargin ~= 3
    error('Incorrect number of arguments')
end

subname = 'grid2fvcom';

global ftbverbose;
if(ftbverbose)
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end

%--------------------------------------------------------------------------
% Get the relevant bits from the FVCOM mesh object
%--------------------------------------------------------------------------
x   = Mobj.x;
y   = Mobj.y;
nVerts = Mobj.nVerts;
nElems = Mobj.nElems;
if(ftbverbose);
    fprintf('info for FVCOM domain\n');
    fprintf('number of nodes: %d\n',nVerts);
    fprintf('number of elems: %d\n',nElems);
end

xc = nodes2elems(x, Mobj);
yc = nodes2elems(y, Mobj);

try
    ntimes = numel(wind.time);
catch
    ntimes = numel(wind.(vars{1}).time);
end

% Interpolate supplied regularly gridded data to FVCOM mesh. Use
% TriScatteredInterp to do the interpolation instead of griddata (should be
% faster).
for vv=1:length(vars)
    if strcmpi(vars{vv}, 'lat') || strcmpi(vars{vv}, 'lon') || strcmpi(vars{vv}, 'x') || strcmpi(vars{vv}, 'y') || strcmpi(vars{vv}, 'time')
        fprintf('skipping variable %s\n', vars{vv})
        continue
    else
        % Preallocate the output arrays
        fvcom.(vars{vv}).data = zeros(nElems,ntimes);
        fvcom.(vars{vv}).node = zeros(nVerts,ntimes);

        for i=1:ntimes
            fprintf('interpolating %s, frame %d of %d\n', vars{vv}, i, ntimes);
            currvar =  wind.(vars{1}).data(:, :, 1);
            % griddata way (cubic interpolation)
            %fvcom.(vars{vv}).node(:,i) = griddata(wind.x,wind.y,currvar,x,y,'cubic');
            %fvcom.(vars{vv}).data(:,i) = griddata(wind.x,wind.y,currvar,xc,yc,'cubic');
            % TriScatteredInterp way (with natural neighbour interpolation)
            ftsin = TriScatteredInterp(wind.x(:), wind.y(:), currvar(:), 'natural');
            fvcom.(vars{vv}).node(:,i) = ftsin(x,y);
            fvcom.(vars{vv}).data(:,i) = ftsin(xc,yc);
        end
        fprintf('interpolation of %s complete\n', vars{vv});
    end
end

