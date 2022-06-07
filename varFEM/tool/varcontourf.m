function [c,h] = varcontourf(x,y,node,elem,uh,n)

if nargin==5, n = 20; end

[x,y] = meshgrid(x,y);
uxy = varInterpolant2d(x,y,node,elem,uh);
[c,h] = contourf(x,y,uxy,n); axis equal
color = jet(10);  colormap(color);