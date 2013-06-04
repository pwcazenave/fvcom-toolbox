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
%
%==========================================================================

matlabrc
close all
clc

addpath('/users/modellers/pica/Code/fvcom-toolbox/utilities')
addpath('/users/modellers/pica/Code/fvcom-toolbox/fvcom_prepro/')

% Temporary paths when on Riqui's machine
addpath('/tmp/pica/fvcom-toolbox/fvcom_prepro/')
addpath('/tmp/pica/fvcom-toolbox/utilities/')

load('/tmp/pica/fvcom-toolbox/tests/data/get_POLCOMS_tsobc_data.mat');

% Perform the interpolation using the new routine.
obc_ts = {'/tmp/pica/fvcom-toolbox/tests/data/Daily.PolcomsErsem.2001.01.nc'};
Mobj_new = get_POLCOMS_tsobc(Mobj, obc_ts);

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
            [~, results.(fnames{ff}).origNodeTimes] = ...
                size(Mobj.(fnames{ff}));
            [results.(fnames{ff}).nNodes, ...
                results.(fnames{ff}).nNodeTimes] = ...
                size(Mobj_new.(fnames{ff}));

            if results.(fnames{ff}).nNodes == Mobj.nVerts
                results.(fnames{ff}).nodeNumber = 'PASS';
            end
            if results.(fnames{ff}).nNodeTimes == ...
                    results.(fnames{ff}).origNodeTimes
                results.(fnames{ff}).numNodeTimes = 'PASS';
            end

            %--------------------------------------------------------------
            % Check the values in the node and element arrays match to
            % reference values.
            %--------------------------------------------------------------
            results.(fnames{ff}).nodeDiff = ...
                Mobj.(fnames{ff}) - ...
                Mobj_new.(fnames{ff});

            results.(fnames{ff}).nodeRange = ...
                max(results.(fnames{ff}).nodeDiff(:)) - ...
                min(results.(fnames{ff}).nodeDiff(:));

            if nodeRange == 0
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
                        results.(fnames{ff}).nNodes)
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
                        results.(fnames{ff}).origNodeTimes, ...
                        results.(fnames{ff}).nNodeTimes)
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
