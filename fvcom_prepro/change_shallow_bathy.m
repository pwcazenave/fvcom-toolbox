function M  = change_shallow_bathy(M, min_depth)
% Deepens shallow nodes by setting a minimum depth and making nodes that
% are too shallow the mean of the surrouding deep nodes. Loop until all
% nodes are deeper than the minuimum depth required.
%
% function [M]=change_shallow_bathy(M, min_depth)
%
% DESCRIPTION:
%
% INPUT:
%  M         = Mesh object
%  min_depth = the minimum depth of nodes
%
% OUTPUT:
%  Nested  Mesh object with altered bathymetry.
%
% EXAMPLE USAGE:
%
%
% Author(s):
%   Rory O'Hara Murray (Marine Scotland Science)
%
% Revision history:
% 2014 sometime - first version
%
%==========================================================================

count = 0;
h=-99; % just to start

while sum(h==-99)>0
    count=count+1;
    disp(['iteration # ' num2str(count)])
    I = M.h<=min_depth;         % find shallow nodes as depth is +ve down.
    If = find(I);
    h = [];
    
    % loop through all shallow nodes with depth<=min_depth
    for ii=1:length(If)
        
        % find elements surrounding the shallow node
        test = [];
        for jj=1:3
            test = [test; find(M.tri(:,jj)==If(ii))];
        end
        % find nodes for all these elements surrounding the shallow node
        % (the surrounding nodes)
        nodes = unique(M.tri(test,:));
        
        htmp = M.h(nodes);
        bla = htmp>min_depth;       % find the nodes that are deeper than min_depth
        if sum(bla)>0
            h(ii) = mean(htmp(bla));% make the shallow node the mean of the surrounding nodes deeper than min_depth
        else
            h(ii) = -99;  % id no deep surroundign nodes then make it -99 and try again.
        end
    end
    M.h(I) = h;  % save changes to the Mobj array.
end

end