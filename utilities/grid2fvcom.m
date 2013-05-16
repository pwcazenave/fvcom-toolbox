function fvcom = grid2fvcom(Mobj,vars,data)
% Interpolate regularly gridded surface forcing data onto a given FVCOM
% grid.
%
% grid2fvcom(Mobj,vars,data)
%
% DESCRIPTION:
%   Takes a given NCEP reanalysis grid file and interpolates the U10 and
%   V10 values onto the specified FVCOM grid file.
%
% INPUT:
%   Mobj - MATLAB mesh object with the following fields:
%       x, y, lon, lat - cartesian and spherical node coordinates. These
%       are transferred to the NetCDF file only and are not used in the
%       interpolation at all.
%       nVerts - number of vertices (nodes) in the unstructured grid.
%       nElems - number of elements in the unstructured grid.
%   vars - a cell array of the variable names to be interpolated on the
%       FVCOM grid in Mobj (e.g. uwnd, U10, vwnd, V10 etc.).
%   data - a struct which contains the following arrays:
%       x - x data (probably best in cartesian for the interpolation)
%       y - y data (probably best in cartesian for the interpolation)
%       The struct must also contain all the variables defined in vars.
%       time - time vector (in Modified Julian Days). If you're using some
%       of the NCEP surface products (e.g. relative humitidy, sea level
%       pressure), you need to supply x and y coordinates for their grids
%       as .xalt and .yalt).
%
% OUTPUT:
%   fvcom - struct of the interpolated data values at the model nodes and
%       element centres. Also includes any variables which were in the
%       input struct but which have not been interpolated (e.g. time).
%
% EXAMPLE USAGE:
%   interpfields = {'uwnd', 'vwnd', 'slp', 'nshf', 'nlwrs', 'nswrs', ...
%       'P_E', 'Et', 'time', 'lon', 'lat', 'x', 'y'};
%   forcing_interp = grid2fvcom(Mobj, interpfields, forcing);
%
% NOTE:
%   The shape of the returned arrays for rhum and slp (via
%   get_NCEP_forcing.m) have sometimes differed from the other vairables
%   (they appear to be projected onto a different grid). Unless you
%   desperately need them, I would suggest omitting them from the
%   interpolation here as this assumes the arrays are all the same size.
%   Alternatively, give data.xalt and data.yalt to specify the alternative
%   grid.
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
%   2013-02-14 Add support for interpolating data on two different grids
%   through the .xalt and .yalt fields in the input data structure. This is
%   handy if you've got data from NCEP from both the Surface and Gaussian
%   products, each of which are on different grids.
%   2013-05-16 Add parallel for loops if the Parallel Computing Toolbox is
%   available (MATLAB parfor loops fail gracefully to for loops if it is
%   not available, in which case no harm, no foul).
%
%==========================================================================

if nargin ~= 3
    error('Incorrect number of arguments')
end

subname = 'grid2fvcom';

global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

% Run jobs on multiple workers if we have that functionality. Not sure if
% it's necessary, but check we have the Parallel Toolbox first.
wasOpened = false;
if license('test', 'Distrib_Computing_Toolbox')
    % We have the Parallel Computing Toolbox, so launch a bunch of workers.
    if matlabpool('size') == 0
        % Force pool to be local in case we have remote pools available.
        matlabpool open local
        wasOpened = true;
    end
end

%--------------------------------------------------------------------------
% Get the relevant bits from the FVCOM mesh object
%--------------------------------------------------------------------------
x = Mobj.x;
y = Mobj.y;
nVerts = Mobj.nVerts;
nElems = Mobj.nElems;
if ftbverbose
    fprintf('info for FVCOM domain\n');
    fprintf('number of nodes: %d\n', nVerts);
    fprintf('number of elems: %d\n', nElems);
end

xc = nodes2elems(x, Mobj);
yc = nodes2elems(y, Mobj);

try
    ntimes = numel(data.time);
catch
    ntimes = numel(data.(vars{1}).time);
end

% Interpolate supplied regularly gridded data to FVCOM mesh. Use
% TriScatteredInterp to do the interpolation instead of griddata (should be
% faster).
for vv = 1:length(vars)
    if strcmpi(vars{vv}, 'time')
        fprintf('transferring variable %s as is\n', vars{vv})
        fvcom.(vars{vv}) = data.(vars{vv});
        continue
    elseif strcmpi(vars{vv}, 'lat') || strcmpi(vars{vv}, 'lon') || strcmpi(vars{vv}, 'x') || strcmpi(vars{vv}, 'y')
        fprintf('reassigning variable %s from unstructured grid\n', vars{vv})
        fvcom.(vars{vv}) = Mobj.(vars{vv});
    else
        % Preallocate the output arrays.
        % Serial version:
        % fvcom.(vars{vv}).data = zeros(nElems, ntimes);
        % fvcom.(vars{vv}).node = zeros(nVerts, ntimes);
        % Also create temporary arrays for the inner loop to be
        % parallelisable (is that a word?):
        tmp_fvcom_data = zeros(nElems, ntimes);
        tmp_fvcom_node = zeros(nVerts, ntimes);
        tmp_data_data = data.(vars{vv}).data; % input to the interpolation

        % Use a parallel loop for the number of time steps we're
        % interpolating (should be quicker, but will use more memory...).
        parfor i = 1:ntimes
            fprintf('interpolating %s, frame %d of %d\n', vars{vv}, i, ntimes);

            % Serial version:
            % currvar = data.(vars{vv}).data(:, :, i);
            % Parallel version:
            currvar = tmp_data_data(:, :, i);

            % griddata way (cubic interpolation)
            %fvcom.(vars{vv}).node(:,i) = griddata(wind.x,wind.y,currvar,x,y,'cubic');
            %fvcom.(vars{vv}).data(:,i) = griddata(wind.x,wind.y,currvar,xc,yc,'cubic');

            % TriScatteredInterp way (with natural neighbour interpolation)
            % Use a try/catch to account for the different grids over which
            % the humidity and sealevel pressure data are sampled.
            try
                ftsin = TriScatteredInterp(data.x(:), data.y(:), currvar(:), 'natural');
            catch err
                warning([err.identifier, ': Some NCEP data are projected' ...
                    ' onto a different grid. Check you have specified' ...
                    ' data.xalt and data.yalt arrays which are on the' ...
                    ' same grid as the data to be interpolated.'])
                ftsin = TriScatteredInterp(data.xalt(:), data.yalt(:), currvar(:), 'natural');
            end
            % Serial version:
            % fvcom.(vars{vv}).node(:,i) = ftsin(x,y);
            % fvcom.(vars{vv}).data(:,i) = ftsin(xc,yc);
            % nnans(1) = sum(isnan(fvcom.(vars{vv}).node(:,i)));
            % nnans(2) = sum(isnan(fvcom.(vars{vv}).data(:,i)));
            % Parallel version:
            tmp_fvcom_node(:, i) = ftsin(x,y);
            tmp_fvcom_data(:, i) = ftsin(xc,yc);
            nnans1 = sum(isnan(tmp_fvcom_node(:, i)));
            nnans2 = sum(isnan(tmp_fvcom_data(:, i)));
            if  nnans1 > 0
                warning('%i NaNs in the interpolated node data. This won''t work with FVCOM.', nnans(1))
            end
            if nnans2 > 0
                warning('%i NaNs in the interpolated element data. This won''t work with FVCOM.', nnans(2))
            end
        end
        % Transfer the temporary arrays back to the relevant struct and
        % clear out the temporary arrays.
        fvcom.(vars{vv}).node = tmp_fvcom_node;
        fvcom.(vars{vv}).data = tmp_fvcom_data;
        clear nnans* tmp_*

        fprintf('interpolation of %s complete\n', vars{vv});
    end
end

if wasOpened
    matlabpool close
end

if ftbverbose
    fprintf('end   : %s \n', subname)
end
