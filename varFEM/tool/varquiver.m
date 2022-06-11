function varquiver(x,y,node,elem,uh,varargin)
% uh = [u, v]

[x,y] = meshgrid(x,y);
u = pdeInterpolant2d(x,y,node,elem,uh(:,1));
v = pdeInterpolant2d(x,y,node,elem,uh(:,2));
quiver(x(:),y(:),u(:),v(:),varargin{:});