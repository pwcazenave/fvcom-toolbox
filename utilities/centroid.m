function [x0,y0,a] = centroid(x,y)
% x and y must be non-scalar vectors of equal sizes defining a polygonal
% region in counterclockwise sequence. a is returned with the polygon's
% area, and (x0,y0) with the coordinates of the polygon's centroid. Note:
% with clockwise order, x0 and y0 are still correct but a is the negative
% of the area.
% 
% RAS - 1/17/05

[m1,n1] = size(x); [m2,n2] = size(y);
n = max(m1,n1);
if [m1,n1] ~= [m2,n2] || min(m1,n1) ~= 1 || n <= 1
    error('Args must be equal-sized non-scalar vectors')
end
x = x(:); y = y(:);
x2 = [x(2:n);x(1)];
y2 = [y(2:n);y(1)];
a = 1/2*sum (x.*y2-x2.*y);
x0 = 1/6*sum((x.*y2-x2.*y).*(x+x2))/a;
y0 = 1/6*sum((x.*y2-x2.*y).*(y+y2))/a;