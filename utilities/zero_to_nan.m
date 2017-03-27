function array = zero_to_nan(array)
%
% Replaces values of zero with NaN.
%
% Author(s)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2017-03-27 Removed the loops and used logical indexing instead.

array(array == 0) = nan;
