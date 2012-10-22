function year = get_NCEP_year(file)
% Extract the year from a give NCEP file name (either 'uwnd.sig995.YYYY.nc'
% or 'vwnd.sig995.YYYY.nc'). Files are those downloaded from
% ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface/.
%
% INPUT:
%   NCEP NetCDF filename (and path)
%
% OUTPUT:
%   NCEP data year
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history:
%   2012-10-15 First version
%
%==========================================================================

[~, tmp_year, ~] = fileparts(file);
tmp_year = regexp(tmp_year, '\.', 'split');
year = str2double(tmp_year(end));
