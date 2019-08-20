function u = Poisson(node,elem,bdFlag,pde)
% Poisson: solve Poisson equation with P1 linear element.
%
%    - Delta u = f   in \Omega, with
%    Dirichlet boundary conditions u = g_D on \Gamma_D,
%    Neumann boundary conditions   grad(u)*n = g_N on \Gamma_N.

% ------------------------- PDE data ------------------------
% pde = struct('uexact',@uexact, 'f',@f, 'g_N',@g_N,  'g_D',@g_D);
% pde = pdedata();

eD = bdFlag.eD;  elemN = bdFlag.elemN;

% -------------- Stiff matrix and load vector -------------
N = size(node,1);  NT = size(elem,1); Ndof = 3;

% area of all triangules
x1 = node(elem(:,1),1); y1 = node(elem(:,1),2);
x2 = node(elem(:,2),1); y2 = node(elem(:,2),2);
x3 = node(elem(:,3),1); y3 = node(elem(:,3),2);
A = 0.5*((x2.*y3-x3.*y2)-(x1.*y3-x3.*y1)+(x1.*y2-x2.*y1));

% All element matrices
k11 = (x3-x2).^2+(y2-y3).^2;
k12 = (x3-x2).*(x1-x3)+(y2-y3).*(y3-y1);
k13 = (x3-x2).*(x2-x1)+(y2-y3).*(y1-y2);
k21 = k12;
k22 = (x1-x3).^2+(y3-y1).^2;
k23 = (x1-x3).*(x2-x1)+(y3-y1).*(y1-y2);
k31 = k13;
k32 = k23;
k33 = (x2-x1).^2+(y1-y2).^2;
K = [k11,k12,k13,k21,k22,k23,k31,k32,k33];
K = K./(4*repmat(A,1,Ndof^2));

% All vectors
f = pde.f;
xc = 1/3*(x1+x2+x3); yc = 1/3*(y1+y2+y3); pc = [xc,yc];
F1 = f(pc).*A./3; F2 = F1; F3 = F1;
F = [F1,F2,F3];

% local --> global
nnz = NT*Ndof^2;
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = K(:);
id = 0;
for i = 1:Ndof
    for j = 1:Ndof
        ii(id+1:id+NT) = elem(:,i);   % zi
        jj(id+1:id+NT) = elem(:,j);   % zj
        id = id + NT;
    end
end

% stiff matrix and load vector
kk = sparse(ii,jj,ss,N,N);
ff = accumarray(elem(:), F(:),[N 1]);

% ------- Neumann boundary conditions -----------
if ~isempty(elemN)
    z1 = node(elemN(:,1),:); % starting points
    z2 = node(elemN(:,2),:); % ending points
    h = sqrt(sum((z2-z1).^2,2));
    g_N = pde.g_N;
    F1 = h./2.*g_N(z1); F2 = h./2.*g_N(z2);
    FN = [F1,F2];
    ff = ff + accumarray(elemN(:), FN(:),[N 1]);
end

% --------- Dirichlet boundary conditions ---------------
if ~isempty(eD)
    g_D = pde.g_D;
    N = size(node,1);
    isBdNode = false(N,1); isBdNode(eD) = true;
    bdNode = find(isBdNode); freeNode = find(~isBdNode);
    pD = node(bdNode,:);
    u = zeros(N,1); u(bdNode) = g_D(pD);
    ff = ff - kk*u;
end

% ------------ Solve the linear system -----------
u(freeNode) = kk(freeNode,freeNode)\ff(freeNode);