function era = get_ERA_forcing(Mobj, files, varargin)
% Reads in ERA-20C netCDF files and outputs a struct containing the
% requested variables for a given spatial domain and time period.
%
% ERA = get_ERA_forcing(Mobj, files, 'varlist', {'lon', 'lat', 'ssr'})
%
% DESCRIPTION:
%   For the given netCDF(s), load the data into a struct for the variables
%   in varlist. If varlist is omitted, all variables are loaded.
%
% INPUT:
%   Mobj - Mesh object containing the following fields:
%       * nativeCoords - 'spherical' or 'cartesian'
%       * lon, lat or x, y - node coordinates (depending on nativeCoords)
%   files - array of paths to the ERA-20C netCDF files
%   Optional keyword/value pairs:
%   'varlist' - give an array of variables to extract from the netCDF
%       files. Missing variables will cause an error.
%   'clip' - set to true for clipping the data to the model domain, false
%       to load the full data in the netCDFs.
%
% OUTPUT:
%   era - struct containing the variables specified in varlist.
%
% NOTES:
%   This script has been written against the netCDFs of the ECMWF-ERA20C
%   data created by the python script ecmwf-era20c.py.
%
% Author(s)
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2015-08-13 First version based loosely on read_ERA_wind.m in the
%   fvcom-toolbox.
%
%==========================================================================

subname = 'get_ERA_forcing';

global ftbverbose
if ftbverbose
	fprintf('\nbegin : %s \n', subname)
end

assert(nargin >= 2, 'Incorrect number of input arguments (2 minimum).')

varlist = [];
clip = false;
for var = 1:2:length(varargin)
    switch varargin{var}
        case 'varlist'
            varlist = varargin{var + 1};
        case 'clip'
            clip = varargin{var + 1};
    end
end

nf = length(files);

for f = 1:nf
    
    assert(exist(files{f}, 'file') == 2, 'file: %s does not exist', files{f})
    
    if ftbverbose
        fprintf('Loading file %s (%i of %i)\n', files{f}, f, nf)
    end

    if isempty(varlist)
        ncdetails = ncinfo(files{f});
        varlist = {ncdetails.Variables.Name};
    end

    % Only do the spatial data on the first file (the data are global).
    if f == 1
        lon = double(ncread(files{f}, 'lon'));
        lat = double(ncread(files{f}, 'lat'));
        lonvector = unique(lon);
        latvector = unique(lat);
        dx = unique(diff(lonvector));
        dy = mean(unique(diff(latvector)));

        if clip
            % To save on memory, cull the data which doesn't come from the
            % domain we've specified (with a 2 element buffer).
            if min(Mobj.lon) < 0 && max(Mobj.lon) < 0
                % Easy, just shift by 360.
                idx_lon = find(lonvector > (min(Mobj.lon) - (2 * dx)) + 360 ...
                    & lonvector < (max(Mobj.lon) + (2 * dx)) + 360);
            elseif min(Mobj.lon) < 0 && max(Mobj.lon) > 0
                % This is the tricky one. We'll do two passes to extract
                % the western chunk first (min(Mobj.lon) + 360 to 360),
                % then the eastern chunk (0 - max(Mobj.lon))
                idx_lon{1} = find(lonvector >= (min(Mobj.lon) - (2 * dx)) + 360)';
                idx_lon{2} = find(lonvector < (max(Mobj.lon) + (2 * dx)))';
                idx_lon = [idx_lon{:}];
            else
                % Dead easy, we're in the eastern hemisphere, so nothing
                % too strenuous here.
                idx_lon = find(lonvector > (min(Mobj.lon) - (2 * dx)) ...
                    & lonvector < (max(Mobj.lon) + (2 * dx)));
            end
            % Latitudes are easy because there's only one system.
            idx_lat = find(latvector > (min(Mobj.lat) - (2 * dy)) & latvector < (max(Mobj.lat) + (2 * dy)));

        else
            idx_lon = 1:length(lonvector);
            idx_lat = 1:length(latvector);
        end
        % Make a grid of the domain
        [era.lon, era.lat] = deal(lonvector(idx_lon), latvector(idx_lat));
        era.lon(era.lon > 180) = era.lon(era.lon > 180) - 360;
    end

    if ftbverbose
        fprintf('Variables:\n')
    end
    for v = 1:length(varlist)
        if ftbverbose
            fprintf('\t%s... ', varlist{v})
        end
        vardump = double(ncread(files{f}, varlist{v}));
        era.(varlist{v}).data = vardump(idx_lon, idx_lat, :);
        era.(varlist{v}).lat = latvector(idx_lat);
        era.(varlist{v}).lon = lonvector(idx_lon);
        era.(varlist{v}).time = double(ncread(files{f}, 'time'));
        clear vardump
        if ftbverbose
            fprintf('done.\n')
        end
    end
end

if ftbverbose
    fprintf('end   : %s \n', subname)
end
