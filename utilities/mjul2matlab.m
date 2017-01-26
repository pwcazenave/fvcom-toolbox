function [ sdOut ] = mjul2matlab( mjul )
%MJUL2MATLAB Convert a Modified Julian date to a Matlab Serial Date
%   Trust Mathworks to use their own date standard when perfectly good ones
%   already existed :-) Still, they're at least easy to convert

% Modified Julian epoch starts at [1858 11 17 0 0 0]. So datenum(that)
% gives the value below.

sdOut = mjul + 678942;

end

