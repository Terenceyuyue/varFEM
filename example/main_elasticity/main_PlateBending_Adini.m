
clc;clear;close all
% -------- Mesh and boudary conditions --------
a1 = 0; b1 = .1; a2 = 0; b2 = 1;
Nx = 10; Ny = 10; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2,'rectangle');
bdStruct = setboundary(node,elem);
% figure,showmesh(node,elem);
% findelem(node,elem);
% findedge(node,elem); findnode(node);

% ------------ PDE data ------------
pde = PlateBendingData;
para = pde.para; D = para.D; nu = para.nu;
f = pde.f;

% -------- Sparse assembling indices -----------
N = size(node,1); NT = size(elem,1); Ndof = 4*3;
elem2 = zeros(NT,Ndof);
elem2(:,[1 4 7 10]) = elem;
elem2(:,[2 5 8 11]) = elem + N;
elem2(:,[3 6 9 12]) = elem + 2*N;

nnz = NT*Ndof^2;
ii = zeros(nnz,1); jj = zeros(nnz,1);
id = 0;
for i = 1:Ndof
    for j = 1:Ndof
        ii(id+1:id+NT) = elem2(:,i);
        jj(id+1:id+NT) = elem2(:,j);
        id = id + NT;
    end
end

% ---------------- parameters used in computation -----------------------
xi = [-1 1 1 -1]; eta = [-1 -1 1 1]; % vertices of [-1,1]^2
x1 = node(elem(:,1),1); y1 = node(elem(:,1),2);
x2 = node(elem(:,2),1); y2 = node(elem(:,2),2);
x3 = node(elem(:,3),1); y3 = node(elem(:,3),2);
x4 = node(elem(:,4),1); y4 = node(elem(:,4),2);
aelem = x2-x1;  belem = y4-y1;


% ----------------- Stiffness matrix ----------------
K = zeros(NT,Ndof^2);
for iel = 1:NT
    Ke = cell(4,4); a = aelem(iel); b = belem(iel);
    for i = 1:4
        for j = 1:4
            xis = xi(i)*xi(j); etas = eta(i)*eta(j);
            Kij = zeros(3,3);
            Kij(1,1) = b^2/a^2*xis*(3+etas) + a^2/b^2*etas*(3+xis) + 2*(7-2*nu)/5*xis*etas;
            Kij(1,2) = -b^2/(2*a)*xi(i)*(3+etas)-a/2*etas*(nu*xi(j)+(2+3*nu)/5*xi(i));
            Kij(1,3) = -a^2/(2*b)*eta(i)*(3+xis)-b/2*xis*(nu*eta(j)+(2+3*nu)/5*eta(i));
            Kij(2,1) = -b^2/(2*a)*xi(j)*(3+etas)-a/2*etas*(nu*xi(i)+(2+3*nu)/5*eta(j));
            Kij(2,2) = b^2/12*(3+xis)*(3+etas)+(1-nu)/30*a^2*etas*(5*etas+3);
            Kij(2,3) = a*b/4*(xi(i)+xi(j))*(eta(i)+eta(j));
            Kij(3,1) = -a^2/(2*b)*eta(j)*(3+xis)-b/2*xis*(nu*eta(i)+(2+3*nu)/5*eta(j));
            Kij(3,2) = a*b/4*(xi(i)+xi(j))*(eta(i)+eta(j));
            Kij(3,3) = a^2/12*(3+xis)*(3+etas)+(1-nu)/30*b^2*xis*(5*etas+3);
            Ke{i,j} = Kij;
        end
    end
    Ke = D/(a*b)*cell2mat(Ke);
    K(iel,:) = reshape(Ke',1,[]);
end
kk = sparse(ii,jj,K(:),3*N,3*N);

% ----------------- load vector ------------------
F = zeros(NT,Ndof);
aux = auxgeometry(node,elem); Centroid = aux.Centroid;
for iel = 1:NT
    pc = Centroid(iel,:);
    a = aelem(iel); b = belem(iel);
    F(iel,:) = f(pc)*a*b/4*[1,a/6,b/6, 1,-a/6,b/6, 1,-a/6,-b/6, 1,a/6,-b/6];
end
ff = accumarray(elem2(:), F(:), [3*N 1]);

% ------------ Dirichlet boundary conditions ----------------
eD = bdStruct.eD;
id = [eD; eD+N; eD+2*N];
isBdDof = false(3*N,1); isBdDof(id) = true;
bdDof = find(isBdDof); freeDof = find(~isBdDof);

g_D = pde.g_D;  pD = node(eD,:);  Dw = pde.Dw;
wD = g_D(pD);  Dw = Dw(pD); wxD = Dw(:,1); wyD = Dw(:,2);

w = zeros(3*N,1); w(bdDof) = [wD;wxD;wyD];

ff = ff - kk*w;

% ------------------ Solver -------------------
w(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
w = w(1:N);


% --------- error analysis -------
uexact = pde.uexact;  we = uexact(node);
figure,
subplot(1,2,1), showsolution(node,elem,w);
zlabel('w');
subplot(1,2,2), showsolution(node,elem,we);
zlabel('we');
Eabs = w-we;  % Absolute errors
figure,showsolution(node,elem,Eabs); zlim('auto');

format shorte
Err = norm(Eabs)./norm(we)





