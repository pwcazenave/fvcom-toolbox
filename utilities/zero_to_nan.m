function [array] = zero_to_nan(array)
%
% Replaces data outside error bands with an interpolated value
%
[rr,cc]=size(array);
for ii=1:cc
ix = find((array(:,ii))==0.);
if ~isempty(ix)
    for i = 1:length(ix); 
     array(ix(i),ii) = NaN;
    end
end
end
