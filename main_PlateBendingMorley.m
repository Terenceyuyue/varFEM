clc;clear;close all
% -------- Mesh and boudary conditions --------
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
Nx = 10; Ny = 10; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);
bdStruct = setboundary(node,elem);
%figure,showmesh(node,elem); findedge(node,elem); findnode(node)

% ------------ PDE data ------------
pde = PlateBendingData1;
para = pde.para;
f = pde.f;

% -------- Sparse assembling indices -----------
N = size(node,1); NT = size(elem,1); Ndof = 6;
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; edge = auxT.edge;
NE = size(edge,1);
elem2 = [elem, elem2edge+N];

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

% -------------- Compute second derivatives of phi_i --------------------
auxG = auxgeometry(node,elem);  area = auxG.area;
x1 = node(elem(:,1),1); y1 = node(elem(:,1),2);
x2 = node(elem(:,2),1); y2 = node(elem(:,2),2);
x3 = node(elem(:,3),1); y3 = node(elem(:,3),2);
% xi, eta
xi = [x2-x3, x3-x1, x1-x2]; eta = [y2-y3, y3-y1, y1-y2];
% edge length
z1 = node(edge(:,1),:); z2 = node(edge(:,2),:);
he = sqrt(sum((z2-z1).^2,2)); L = he(elem2edge);
% sij
sb = cell(3,3);
for i = 1:3
    for j = 1:3
        sb{i,j} = -xi(:,i).*xi(:,j) - eta(:,i).*eta(:,j);
    end
end
% second derivatives of phi_i
b11 = zeros(NT,6); b22 = b11; b12 = b11; % [phi_i, i=1:6]
ind = [1 2 3; 2 3 1; 3 1 2];  % cyclic permutation
for i = 1:3
    j = ind(i,2); k = ind(i,3);
    c0 = 1./(2*area.^2);  c1 = sb{i,j}./L(:,j).^2;  c2 = sb{k,i}./L(:,k).^2;
    % i = 1,2,3
    b11(:,i) = c0.* (eta(:,i).^2 - c1.*eta(:,j).^2 - c2.*eta(:,k).^2);
    b22(:,i) = c0.* (xi(:,i).^2 - c1.*xi(:,j).^2 - c2.*xi(:,k).^2);
    b12(:,i) = -c0.* (xi(:,i).*eta(:,i) - c1.*xi(:,j).*eta(:,j) - c2.*xi(:,k).*eta(:,k));
    % i = 4,5,6
    ci = 1./(area.*L(:,i));
    b11(:,3+i) = ci.*eta(:,i).^2;
    b22(:,3+i) = ci.*xi(:,i).^2;
    b12(:,3+i) = -ci.*xi(:,i).*eta(:,i);
end

% ----------- First stiffness matrix -----------
K = zeros(NT,Ndof^2);  s = 1;
for i = 1:Ndof
    for j = 1:Ndof
        K(:,s) = b11(:,i).*b11(:,j) + b22(:,i).*b22(:,j) ...
            + para.nu*(b11(:,i).*b22(:,j) + b22(:,i).*b11(:,j)) ...
            + 2*(1-para.nu)*b12(:,i).*b12(:,j);
        s = s+1;
    end
end
K = para.D*area.*K;

% ----------- Second stiffness matrix -----------
% Gauss quadrature rule
[lambda,weight] = quadpts(3);
G = zeros(NT,Ndof^2); F = zeros(NT,Ndof);
base = zeros(length(lambda),Ndof);
if isa(para.c,'double'), cf = @(xy) para.c+0*xy(:,1); end
for iel = 1:NT
    vK = node(elem(iel,:),:); % vertices of K
    xy = lambda*vK;
    cxy = cf(xy); fxy = f(xy);
    for i = 1:3
        j = ind(i,2); k = ind(i,3);
        c1 = sb{i,j}(iel)/L(iel,j)^2;
        c2 = sb{k,i}(iel)/L(iel,k)^2;   
        c3 = c1+c2;     
        % i = 1,2,3
        base(:,i) = lambda(:,i).^2 + c3*lambda(:,j).*lambda(:,k) ...
            + c2*lambda(:,k).*lambda(:,i) + c1*lambda(:,i).*lambda(:,j);
        % i = 4,5,6
        ci = 2*area(iel)/L(iel,i);
        base(:,3+i) = ci*lambda(:,i).*(lambda(:,i)-1);
    end
    
    s = 1;
    for i = 1:Ndof
        for j = 1:Ndof
            gs = cxy.*base(:,i).*base(:,j);
            G(iel,s) = area(iel)*weight*gs;
            s = s+1;
        end
    end
    
    fv = fxy.*base;
    F(iel,:) = area(iel)*weight*fv;
end

kk = sparse(ii,jj,K(:)+G(:),N+NE,N+NE);

% ------------- Assemble load vector ------------
ff = accumarray(elem2(:), F(:), [N+NE 1]);

% ------------ Dirichlet boundary conditions ----------------
eD = bdStruct.eD; elemD = bdStruct.elemD;
bdIndex = bdStruct.bdIndex; % indices of boundary edges
id = [eD; bdIndex+N];
isBdDof = false(N+NE,1); isBdDof(id) = true;
bdDof = find(isBdDof); freeDof = find(~isBdDof);

g_D = pde.g_D;  
pD = node(eD,:); wD = g_D(pD);


Dw = pde.Dw;
z1 = node(elemD(:,1),:); z2 = node(elemD(:,2),:); zc = (z1+z2)./2;
e = z1-z2;  % e = z2-z1
ne = [-e(:,2),e(:,1)];  ne = ne./he(bdIndex);
wn = sum(Dw(zc).*ne,2);

w = zeros(N+NE,1); w(bdDof) = [wD;wn];
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
zlabel('Err');



