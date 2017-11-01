function Mobj = get_HYCOM_tsobc(Mobj, hycom, varlist)
% Read temperature, salinity, u and v data from the HYCOM model output
% structure and interpolate onto the open boundaries in Mobj.
%
% function Mobj = get_HYCOM_tsobc(Mobj, hycom, varlist)
%
% DESCRIPTION:
%    Interpolate temperature and salinity values onto the FVCOM open
%    boundaries at all sigma levels.
%
% INPUT:
%   Mobj    = MATLAB mesh structure which must contain:
%               - Mobj.siglayz - sigma layer depths for all model nodes.
%               - Mobj.siglayzc - sigma layer depths for all model elements.
%               - Mobj.lon, Mobj.lat - node coordinates (lat/long).
%               - Mobj.lonc, Mobj.latc - element coordinates (lat/long).
%               - Mobj.read_obc_nodes - cell array of open boundary nodes.
%               - Mobj.read_obc_elems - cell array of open boundary
%               elements (only if using velocities - use find_nested_region
%               to get these indices).
%               - Mobj.h - water depths at nodes.
%               - Mobj.tri - triangulation table for the grid (nv in FVCOM
%               terms).
%               - Mobj.nObcNodes - number of nodes in each open boundary.
%   hycom   = Struct with HYCOM data covering the model domain. Unless
%             varlist is specified (see below), all 4D fields will be
%             interpolated onto the open boundaries (1-3D data will be
%             ignored).
%   varlist = [optional] cell array of variable (field) names to use from
%             hycom.
%
% OUTPUT:
%    Mobj = MATLAB structure with new fields whose names match those given
%    in hycom. The fields have sizes (sum(Mobj.nObcNodes), sigma, time).
%    The time dimension is determined based on the number of time steps in
%    hycom. The ts_time variable is just the input file times in Modified
%    Julian Day.
%
% EXAMPLE USAGE
%    hycom = get_HYCOM_forcing(Mobj, [51500, 51531]);
%    Mobj = get_HYCOM_tsobc(Mobj, hycom, {'u', 'v', 'temperature', 'salinity'})
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2013-09-03 First version based on get_POLCOMS_tsobc.m.
%    2014-04-28 Update interp1 function to use pchip instead of csap as the
%    latter will be removed in a future version of MATLAB and the
%    innumerable warnings were doing my nut in. I checked the output using
%    the new interp1 call and it's identical to the old version. Also
%    update the parallel toolbox stuff for the same reason (future
%    removal).
%    2015-05-21 Remove the old parallel processing bits and replace with
%    the current versions.
%    2016-03-15 Add fallback interpolation to inverse distance weighted if
%    the triangular interpolation fails (which can happen if the points
%    supplied are all in a line, for example).
%    2017-01-27 Subset the coarse data (HYCOM, CMEMS etc.). This yields a
%    significant performance improvement (both in memory usage and time).
%    2017-02-16 Further performance improvement by only using the coarse
%    data in the vicinity of the open boundary nodes.
%    2017-10-12 Fix bug in indexing the open boundary positions which may
%    have caused the interpolation to fail as the identified positions
%    would be too far from the open boundary nodes.
%
%==========================================================================

% This is just a wrapper around the more generic interp_coarse_to_obc
% function. This is to maintain backwards compatibility.
Mobj = interp_coarse_to_obc(Mobj, hycom, varlist);
