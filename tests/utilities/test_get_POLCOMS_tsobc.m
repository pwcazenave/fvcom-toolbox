% Unit test for get_POLCOMS_tsobc.
%
% DESCRIPTION:
%   Currently checks against a reference data set for the following:
%       - number of nodes in the output
%       - number of sigma layers in the output
%       - number of time steps in the output
%       - range of values in the node arrays
%
% It uses a simplified POLCOMS NetCDF file from January, 2001 as the base
% input. The mesh object (Mobj) contains the required input for
% get_POLCOMS_tsobc as well as a set of 'known good' results
% (Mobj.temperature, Mobj.salt and Mobj.ts_times) for comparison against
% the new result.
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2013-05-17 First version.
%   2016-06-02 Updated to actually compare the interpolated data rather
%   than just the grid information.
%
%==========================================================================

matlabrc
close all
clc

% Set up our test environment.
[base, subname] = fileparts(mfilename('fullpath'));
addpath(fullfile(base, '../../fvcom_prepro'))

load(fullfile(base, '../data/get_POLCOMS_tsobc_data.mat'));

% Perform the interpolation using the new routine.
obc_ts = {fullfile(base, '../data/Daily.PolcomsErsem.2001.01.nc')};
Mobj_new = get_POLCOMS_tsobc(Mobj, obc_ts);

% Add the necessary information for the checks.
t = delaunayTriangulation(Mobj_new.lon, Mobj_new.lat);
Mobj_new.tri = t.ConnectivityList;
clear t
Mobj.nVerts = length(Mobj.lon);
Mobj.nElems = size(Mobj.tri, 1);
Mobj_new.nVerts = length(Mobj_new.lon);
Mobj_new.nElems = size(Mobj_new.tri, 1);

% Check we have the temperature, salinity and time fields in the new mesh
% object.
fnames = fieldnames(Mobj_new);

for ff = 1:length(fnames)
    switch fnames{ff}
        case {'temperature', 'salt', 'ts_times'}
            assert(isfield(Mobj_new, fnames{ff}), 'Missing field %s', fnames{ff})
    end
end

% Clear out fields which don't exist in both Mobj and Mobj_new.
fnames = intersect(fnames, fieldnames(Mobj));

%%

results = struct();

for ff = 1:length(fnames)

    results.(fnames{ff}) = struct();

    switch fnames{ff}
        case {'siglayz', 'lon', 'lat', 'obc_nodes', 'nObcNodes', 'ts_times'}

            results.(fnames{ff}).vectorValues = 'FAIL';

            results.(fnames{ff}).check = ...
                Mobj.(fnames{ff}) - Mobj_new.(fnames{ff});
            checkDiff = max(results.(fnames{ff}).check) - ...
                min(results.(fnames{ff}).check);
            if checkDiff == 0
                results.(fnames{ff}).vectorValues = 'PASS';
            end

        otherwise

            %--------------------------------------------------------------
            % Set the pass/fail flags for the tests. Assume fail and only
            % change if proven otherwise.
            %--------------------------------------------------------------
            results.(fnames{ff}).nodeNumber = 'FAIL';
            results.(fnames{ff}).elementNumber = 'FAIL';
            results.(fnames{ff}).numNodeTimes = 'FAIL';
            results.(fnames{ff}).nodeValues = 'FAIL';

            %--------------------------------------------------------------
            % Check we have the same number of points and time steps in the
            % new interpolation as in the original.
            %--------------------------------------------------------------

            % Get number of new nodes, times and elements
            results.(fnames{ff}).nVerts = Mobj_new.nVerts;
            results.(fnames{ff}).nElems = Mobj_new.nElems;
            results.(fnames{ff}).nTimes = length(Mobj_new.ts_times);
            results.(fnames{ff}).orignTimes = length(Mobj.ts_times);

            % Compare old and new nodes, elements and times.
            if results.(fnames{ff}).nVerts == Mobj.nVerts
                results.(fnames{ff}).nodeNumber = 'PASS';
            end
            if results.(fnames{ff}).nTimes == length(Mobj.ts_times)
                results.(fnames{ff}).numNodeTimes = 'PASS';
            end
            if results.(fnames{ff}).nElems == Mobj.nElems
                results.(fnames{ff}).elementNumber = 'PASS';
            end

            %--------------------------------------------------------------
            % Check the values in the node arrays match the reference
            % values.
            %--------------------------------------------------------------
            results.(fnames{ff}).nodeDiff = ...
                Mobj.(fnames{ff}) - ...
                Mobj_new.(fnames{ff});

            results.(fnames{ff}).nodeRange = ...
                max(results.(fnames{ff}).nodeDiff(:));

            if results.(fnames{ff}).nodeRange == 0
                results.(fnames{ff}).nodeValues = 'PASS';
            end
    end
end

%%
%--------------------------------------------------------------------------
% Print a summary of the testing
%--------------------------------------------------------------------------
totalTests = 0;
totalPasses = 0;

for ff = 1:length(fnames)
    resultnames = fieldnames(results.(fnames{ff}));
    numRes = length(resultnames);

    for fi = 1:numRes

        % Skip if the field is not a string (we're only interested in
        % pass/fail really.
        if ~ischar(results.(fnames{ff}).(resultnames{fi}))
            continue
        else
            % Increment the number of tests performed.
            totalTests = totalTests + 1;
        end

        % Get the total number of PASSed tests.
        if strcmp(results.(fnames{ff}).(resultnames{fi}), 'PASS')
            totalPasses = totalPasses + 1;
        end

        S = results.(fnames{ff}).(resultnames{fi});

        switch resultnames{fi}
            case 'vectorValues'
                fprintf('%s %s values test\n', S, fnames{ff})
                if strcmp(S, 'FAIL')
                    fprintf('\tmin/max of %s range: %f, %f\n', ...
                        fnames{ff}, ...
                        min(results.(fnames{ff}).check), ...
                        max(results.(fnames{ff}).check))
                end

            case 'nodeNumber'
                fprintf('%s %s node number test\n', S, fnames{ff})
                if strcmp(S, 'FAIL')
                    fprintf('\toriginal/new number of %s nodes: %d, %d\n', ...
                        fnames{ff}, ...
                        Mobj.nVerts, ...
                        results.(fnames{ff}).nVerts)
                end

            case 'elementNumber'
                fprintf('%s %s element number test\n', S, fnames{ff})
                if strcmp(S, 'FAIL')
                    fprintf('\toriginal/new number of %s elements: %d, %d\n', ...
                        fnames{ff}, ...
                        Mobj.nElems, ...
                        results.(fnames{ff}).nElems)
                end

            case 'numNodeTimes'
                fprintf('%s %s node time steps test\n', S, fnames{ff})
                if strcmp(S, 'FAIL')
                    fprintf('\toriginal/new number of %s node times: %d, %d\n', ...
                        fnames{ff}, ...
                        results.(fnames{ff}).orignTimes, ...
                        results.(fnames{ff}).nTimes)
                end

            case 'numElementTimes'
                fprintf('%s %s element time steps test\n', S, fnames{ff})
                if strcmp(S, 'FAIL')
                    fprintf('\toriginal/new number of %s element times: %d, %d\n', ...
                        fnames{ff}, ...
                        results.(fnames{ff}).origElementTimes, ...
                        results.(fnames{ff}).nElementTimes)
                end

            case 'nodeValues'
                fprintf('%s %s node values test\n', S, fnames{ff})
                if strcmp(S, 'FAIL')
                    fprintf('\trange of %s node values: %d\n', ...
                        fnames{ff}, ...
                        results.(fnames{ff}).nodeRange)
                end

            case 'elementValues'
                fprintf('%s %s element values test\n', S, fnames{ff})
                if strcmp(S, 'FAIL')
                    fprintf('\trange of %s element values: %d\n', ...
                        fnames{ff}, ...
                        results.(fnames{ff}).elemRange)
                end
        end
    end
end

fprintf('\n------------------SUMMARY------------------\n')
fprintf('           %d of %d tests passed', totalPasses, totalTests)
fprintf('\n-------------------------------------------\n')
