function fvcom = grid2fvcom(Mobj, vars, data, varargin)
% Interpolate regularly gridded surface forcing data onto a given FVCOM
% grid.
%
% fvcom = grid2fvcom(Mobj,vars,data)
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
%   varargin - keyword/argument pairs:
%       add_elems - true/false - set to true to enable interpolation of
%       field onto both nodes and elements (defaults to both). This is
%       useful if you run into memory issues as this can effectively halve
%       the memory requirements.
%
% OUTPUT:
%   fvcom - struct of the interpolated data values at the model nodes and
%       element centres (unless add_elems is false, in which case, only
%       nodes). Also includes any variables which were in the input struct
%       but which have not been interpolated (e.g. time).
%
% EXAMPLE USAGE:
%   interpfields = {'uwnd', 'vwnd', 'slp', 'nshf', 'nlwrs', 'nswrs', ...
%       'P_E', 'Et', 'time', 'lon', 'lat', 'x', 'y'};
%   forcing_interp = grid2fvcom(Mobj, interpfields, forcing);
%
% NOTE:
%   The shape of the returned arrays for rhum and slp (via
%   get_NCEP_forcing.m) have sometimes differed from the other vairables
%   (they appear to be projected onto a different grid). Give data.xalt and
%   data.yalt to specify the alternative grid.
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
%   2013-07-18 Add more elegant case statement rather than using string
%   comparisons.
%   2013-08-06 Fix fairly significant bug in which the position arrays were
%   transposed relative to the data arrays. The code now checks for the
%   dimensions and warns if they have been flipped to match. There is no
%   checking that the flip has worked because the xalt and yalt arrays
%   complicate things too much for me to figure out today. If you want to
%   implement that functionality, please do so! I also added commented-out
%   figure at the end to check the interpolation has worked properly,
%   should you wish to check manually.
%   2013-12-04 Check for the presence of the input fields being requested
%   in the input struct to avoid finding out that the last field in vars
%   doesn't exist in data. Change the way the alternative coordinate arrays
%   are used to accommodate subtleties in the parallel code in MATLAB.
%   2015-05-20 Update the parallel processing commands.
%   2015-05-22 Add option to disable output on element centres.
%   2016-06-02 Fix the alternative grid handling in the parallel loop.
%   Remove commented out code too.
%
%==========================================================================

if nargin ~= 3 && isempty(varargin)
    error('Incorrect number of arguments')
end

do_elems = false;
for v = 1:2:length(varargin)
    switch varargin{v}
        case 'add_elems'
            if varargin{v + 1}
                do_elems = true;
            end
    end
end

subname = 'grid2fvcom';

global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

% Before we go too far into this, check we have all the fields in the input
% data that are being requested.
for ff = 1:length(vars)
    assert(isfield(data, vars{ff}), ...
        'Missing field %s in the input data struct.', vars{ff})
end

% Run jobs on multiple workers if we have that functionality. Not sure if
% it's necessary, but check we have the Parallel Toolbox first.
if license('test', 'Distrib_Computing_Toolbox')
    % We have the Parallel Computing Toolbox, so launch a bunch of workers.
    if isempty(gcp('nocreate'))
        % Force pool to be local in case we have remote pools available.
        parpool('local');
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
    switch vars{vv}
        case 'time'
            fprintf('transferring variable %s as is\n', vars{vv})
            fvcom.(vars{vv}) = data.(vars{vv});
            continue

        case {'lat', 'lon', 'x', 'y'}
            fprintf('reassigning variable %s from unstructured grid\n', vars{vv})
            fvcom.(vars{vv}) = Mobj.(vars{vv});

        case {'xalt', 'yalt'}
            % Only exist for the interpolation of some data on an
            % alternative grid. 
            fprintf('skipping %s\n', vars{vv})

        otherwise
            % Preallocate the output arrays. Also create temporary arrays
            % for the inner loop to be parallelisable (is that a word?):
            tmp_fvcom_data = zeros(nElems, ntimes);
            tmp_fvcom_node = zeros(nVerts, ntimes);
            try
                tmp_data_data = data.(vars{vv}).data; % input to the interpolation
            catch msg
                fprintf('Trying for alternative data structure. (%s) ', ...
                    msg.message)
                tmp_data_data = data.(vars{vv}); % input to the interpolation
                fprintf('success!\n')
            end

            xx = data.x(:);
            yy = data.y(:);

            % Check the shapes of the input data match those of the
            % position arrays.
            [fvx, fvy] = size(data.x);
            [ncx, ncy, ~] = size(tmp_data_data);

            if isfield(data, 'xalt')
                [fvxalt, fvyalt] = size(data.xalt);
                if (ncx ~= fvx || ncy ~= fvy) || (ncx ~= fvxalt || ncy ~= fvyalt)
                    % Flipping the input array so it hopefully matches the
                    % position arrays.
                    tmp_data_data = permute(tmp_data_data, [2, 1, 3]);
                    warning('Transposed ''%s'' input data to match position array dimensions', vars{vv})
                end
            else
                if (ncx ~= fvx || ncy ~= fvy)
                    % Flipping the input array so it hopefully matches the
                    % position arrays.
                    tmp_data_data = permute(tmp_data_data, [2, 1, 3]);
                    warning('Transposed ''%s'' input data to match position array dimensions', vars{vv})
                end
                % If we have a land mask, mask off the coastal and land
                % points in the coordinates arrays.
                if isfield(data, 'lsm')
                    xx(data.lsm ~= 0) = [];
                    yy(data.lsm ~= 0) = [];
                end
            end

            % Use a parallel loop for the number of time steps we're
            % interpolating.
            varname = vars{vv};
            parfor i = 1:ntimes
                if ftbverbose
                    fprintf('interpolating %s, frame %d of %d\n', varname, i, ntimes);
                end

                currvar = tmp_data_data(:, :, i);

                % TriScatteredInterp way (with natural neighbour
                % interpolation). Instead of the quite crude try/catch that
                % was here, count the number of elements in the coordinate
                % (xx and yy) and data (currvar) arrays: if they differ,
                % try the same thing with the xxalt and xyalt coordinate
                % arrays. If they still differ, then error out. The reason
                % for this different approach is that the parfor sometimes
                % failed for me when using non-NCEP data as the source for
                % the interpolation.
                ndata = numel(currvar(~isnan(currvar)));
                nxx = numel(xx);
                nyy = numel(yy);
                % Get the alternate grid, if we have it.
                if isfield(data, 'xalt')
                    xxalt = data.xalt(:);
                    yyalt = data.yalt(:);
                    if isfield(data, 'lsmalt')
                        xxalt(data.lsmalt ~= 0) = [];
                        yyalt(data.lsmalt ~= 0) = [];
                    end
                end
                assert(nxx == nyy, 'Inconsistent coordinate array sizes.')
                if nxx == ndata
                    ftsin = TriScatteredInterp(...
                        xx, ...
                        yy, ...
                        currvar(~isnan(currvar(:))), ...
                        'natural');
                elseif isfield(data, 'xalt') && numel(xxalt) == ndata
                    ftsin = TriScatteredInterp(...
                        xxalt, ...
                        yyalt, ...
                        currvar(~isnan(currvar(:))), ...
                        'natural');
                else
                    error('Can''t interpolate the data: non-matching coordinate array sizes.')
                end

                tmp_fvcom_node(:, i) = ftsin(x, y);
                nnans1 = sum(isnan(tmp_fvcom_node(:, i)));
                if  nnans1 > 0
                    warning('%i NaNs in the interpolated node data. This won''t work with FVCOM.', nnans1)
                end
                if do_elems
                    tmp_fvcom_data(:, i) = ftsin(xc, yc);
                    nnans2 = sum(isnan(tmp_fvcom_data(:, i)));
                    if nnans2 > 0
                        warning('%i NaNs in the interpolated element data. This won''t work with FVCOM.', nnans2)
                    end
                end
            end
            % Transfer the temporary arrays back to the relevant struct and
            % clear out the temporary arrays.
            fvcom.(vars{vv}).node = tmp_fvcom_node;
            if do_elems
                fvcom.(vars{vv}).data = tmp_fvcom_data;
            end
            clear nnans* tmp_*

            if ftbverbose
                fprintf('interpolation of %s complete\n', vars{vv});
            end
    end
end

cleaner = onCleanup(@() delete(gcp('nocreate')));

if ftbverbose
    fprintf('end   : %s \n', subname)
end

%% Debugging figure to check the interpolation makes sense.

% close all
%
% figure
%
% vartoplot = 'nswrf';
% tidx = 12; % time index
%
% subplot(2, 1, 1)
% try
%     pcolor(data.lon, data.lat, data.(vartoplot).data(:, :, tidx)')
% catch
%     pcolor(data.lon, data.lat, data.(vartoplot)(:, :, tidx)')
% end
% shading flat
% axis('equal', 'tight')
% title([vartoplot, ' (regularly gridded)'])
% colorbar
% try
%     caxis([min(fvcom.(vartoplot).data(:, tidx)), max(fvcom.(vartoplot).data(:, tidx))])
% catch
%     caxis([min(fvcom.(vartoplot)(:, tidx)), max(fvcom.(vartoplot)(:, tidx))])
% end
% axis([min(Mobj.lon), max(Mobj.lon), min(Mobj.lat), max(Mobj.lat)])
%
% subplot(2, 1, 2)
% try
%     patch('Vertices', [Mobj.lon, Mobj.lat], 'Faces', Mobj.tri, 'cData', fvcom.(vartoplot).data(:, tidx));
% catch
%     patch('Vertices', [Mobj.lon, Mobj.lat], 'Faces', Mobj.tri, 'cData', fvcom.(vartoplot)(:, tidx));
% end
% shading flat
% axis('equal')
% axis([min(data.lon(:)), max(data.lon(:)), min(data.lat(:)), max(data.lat(:))])
% title([vartoplot, ' (interpolated)'])
% colorbar
% try
%     caxis([min(fvcom.(vartoplot).data(:, tidx)), max(fvcom.(vartoplot).data(:, tidx))])
% catch
%     caxis([min(fvcom.(vartoplot)(:, tidx)), max(fvcom.(vartoplot)(:, tidx))])
% end
% axis([min(Mobj.lon), max(Mobj.lon), min(Mobj.lat), max(Mobj.lat)])

