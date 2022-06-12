function uxy = pdeInterpolant2d(x,y,node,elem,u_node)

p = node';
t = elem';
t(4,:) = 1; % the subdomain number
pdeInt = pdeInterpolant(p, t, u_node);

% Compute the interpolated value at coordinates within the mesh

[nx,ny] = size(x);
uxy = zeros(nx,ny);
% for j = 1:ny
%     uxy(:,j) = pdeInt.evaluate(x(:,j), y(:,j));
% end
uxy(:) = pdeInt.evaluate(x(:),y(:));