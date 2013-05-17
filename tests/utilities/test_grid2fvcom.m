% Unit test for grid2fvcom.

% Check we always produce the same interpolated values.

% Use some NCEP data for the Irish Sea as the base input. This data is from
% around January, 2001. This includes an unstructured grid object (Mobj),
% the NCEP forcing data struct (forcing) and a 'known good' result
% (forcing_interp) for comparison against the new result.
load('/tmp/pica/fvcom-toolbox/tests/mat/grid2fvcom_data.mat');

interpfields = {'uwnd', 'vwnd', 'slp', 'nshf', 'nlwrs', 'nswrs', 'P_E', ...
    'Et', 'time', 'lon', 'lat', 'x', 'y'};

% Perform the interpolation using the new routine and check the outputs are
% the same.
forcing_interp_new = grid2fvcom(Mobj, interpfields, forcing);

% Now we have the new results, compare the ranges of values in the old and
% the new to make sure they're the same.
fnames = fieldnames(forcing_interp);
for ff = 1:length(fnames)
    node_diff = forcing_interp.(fnames{ff}).node - ...
        forcing_interp_new.(fnames{ff}).node;
    elem_diff = forcing_interp.(fnames{ff}).data - ...
        forcing_interp_new.(fnames{ff}).data;

    fprintf('%s min/max of node differences: %f\n', min(node_diff), max(node_diff))
    fprintf('%s min/max of element differences: %f\n', min(elem_diff), max(elem_diff))

end
