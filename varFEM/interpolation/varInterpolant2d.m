function uxy = varInterpolant2d(x,y,node,elem,uh)

z1 = node(elem(:,1),:);  x1 = z1(:,1);  y1 = z1(:,2);
z2 = node(elem(:,2),:);  x2 = z2(:,1);  y2 = z2(:,2);
z3 = node(elem(:,3),:);  x3 = z3(:,1);  y3 = z3(:,2);

xi = [x2-x3, x3-x1, x1-x2];
eta = [y2-y3, y3-y1, y1-y2];
omega = [x2.*y3-x3.*y2, x3.*y1-x1.*y3, x1.*y2-x2.*y1];
S = sum(omega,2)/2;

lam1 = @(x,y) 1./(2*S).*(eta(:,1)*x - xi(:,1)*y + omega(:,1));
lam2 = @(x,y) 1./(2*S).*(eta(:,2)*x - xi(:,2)*y + omega(:,2));
lam3 = @(x,y) 1./(2*S).*(eta(:,3)*x - xi(:,3)*y + omega(:,3));

[nx,ny] = size(x);
uxy = NaN(nx,ny);
for i = 1:nx
    for j = 1:ny
        xc = x(i,j); yc = y(i,j);
        phi1 = lam1(xc,yc); 
        phi2 = lam2(xc,yc); 
        phi3 = lam3(xc,yc);
        id = find((phi1>=0) & (phi2>=0) & (phi3>=0));
        if isempty(id), continue; end
        iel = id(1);  
        uxy(i,j) = uh(elem(iel,1))*phi1(iel) ...
            + uh(elem(iel,2))*phi2(iel) ...
            + uh(elem(iel,3))*phi3(iel);
    end
end