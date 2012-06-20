function res = isintriangle(xt,yt,x0,y0)
%  determine if point (x0,y0) is in triangle defined by nodes (xt(3),yt(3))    
%
%  function res = isintriangle(xt,yt,x0,y0)
%
%  determine if point (x0,y0) is in triangle defined by nodes (xt(3),yt(3))    |
%  using algorithm used for scene rendering in computer graphics               |
%  algorithm works well unless particle happens to lie in a line parallel      |
%  to the edge of a triangle.                                                  |
%  This can cause problems if you use a regular grid, say for idealized        |
%  modelling and you happen to seed particles right on edges or parallel to    |
%  edges.                                                                      |
%==============================================================================|

 res = 0;
 f1 = (y0-yt(1))*(xt(2)-xt(1)) - (x0-xt(1))*(yt(2)-yt(1));
 f2 = (y0-yt(3))*(xt(1)-xt(3)) - (x0-xt(3))*(yt(1)-yt(3));
 f3 = (y0-yt(2))*(xt(3)-xt(2)) - (x0-xt(2))*(yt(3)-yt(2));
 if(f1*f3 >= 0.0 & f3*f2 >= 0.0) 
   res = 1;
 end;

 return

