function rh = qair2rh(qair, temp, varargin)
% Convert specific humitidy to relative humidity.
%
% INPUTS:
%   qair - Specific humidity. Dimensionless ratio of water mass /
%       total air mass.
%   temp - Temperature. Celsius.
%   press - Pressure, optional. Millibar. Assume as 1013.25 if
%       missing.
%
% OUTPUT:
%   rh  - relative humidity, ratio of actual water mixing ratio to
%       saturation mixing ratio.
%
% NOTES:
%   Translated from the R code here:
%       http://earthscience.stackexchange.com/questions/2360
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)

if nargin == 2
    press = 1013.25;
end
es = 6.112 * exp((17.67 * temp) / (temp + 243.5));
e = qair * press / (0.378 * qair + 0.622);
rh = e / es;
rh(rh > 1) = 1;
rh(rh < 0) = 0;
