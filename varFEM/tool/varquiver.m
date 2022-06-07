function varquiver(x,y,node,elem,uh,varargin)
% uh = [u, v]

[x,y] = meshgrid(x,y);
u = varInterpolant2d(x,y,node,elem,uh(:,1));
v = varInterpolant2d(x,y,node,elem,uh(:,2));
quiver(x(:),y(:),u(:),v(:),varargin{:});