function [velMin,velMax,velMed,velMean,velStd]=get_velocity(u,v)
% Get velocity stats from u and v components.
% 
% 

vel = sqrt(u(:).^2+v(:).^2);
velMin = min(vel);
velMax = max(vel);
velMed = median(vel);
velMean = mean(vel);
velStd = std(vel);
