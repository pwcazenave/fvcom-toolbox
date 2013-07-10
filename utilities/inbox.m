function [pts] = inbox(box,x,y)
% determine if points lie in a box

pts = find( x > box(1) & x < box(2) & y >box(3) & y < box(4));


