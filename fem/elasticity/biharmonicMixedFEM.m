function [u,w] = biharmonicMixedFEM(node,elem,pde,bdStruct)
% BiharmonicMixedFEM solves the biharmonic equation
%
%     Laplace^2 u = f;   [a1,b1] * [a2,b2]
%               u = g_D;
%            Dn u = g_N.
%
% by writing in a mixed form
%
%       - Laplace u = w
%       - Laplace w = f
%
% Unlike the conforming or nonconforming FEMs, the second boundary condition
% is a Neumann boundary condition in this case.

N = size(node,1); NT = size(elem,1); Ndof = 3;
% -------------- sparse assembling index ----------------
nnz = NT*Ndof^2;
ii = zeros(nnz,1); jj = zeros(nnz,1);
id = 0;
for i = 1:Ndof
    for j = 1:Ndof
        ii(id+1:id+NT) = elem(:,i);   % zi
        jj(id+1:id+NT) = elem(:,j);   % zj
        id = id + NT;
    end
end

ii11 = ii;   jj11 = jj;  ii12 = ii;   jj12 = jj+N;
ii21 = ii+N; jj21 = jj;

% -------------- Stiffness matrix B -------------
% B: (Dphi_i, Dphi_j)_K
[Dphi,area] = gradbasis(node,elem);
B = zeros(NT,Ndof^2); s = 1;
for i = 1:Ndof
    for j = 1:Ndof;
        B(:,s) = sum(Dphi(:,:,i).*Dphi(:,:,j),2).*area;
        s = s+1;
    end
end

% -------------- Stiffness matrix A and load vector -------------
% Gauss quadrature rule
[lambda,weight] = quadpts(3);  nI = length(weight);
% A: -(phi_i, phi_j)_K
A = zeros(NT,Ndof^2); F = zeros(NT,Ndof);
for p = 1:nI
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    % Second stiffness matrix
    s = 1;
    for i = 1:Ndof
        for j = 1:Ndof
            gs = lambda(p,i).*lambda(p,j);
            A(:,s) = A(:,s) + weight(p)*gs;
            s = s+1;
        end
    end
    % load vector
    F = F + weight(p)*pde.f(pxy)*lambda(p,:);
end
A = -repmat(area,1,Ndof^2).*A; F = repmat(area,1,Ndof).*F;

% ------------ Assemble stiffness matrix and load vector ---------
ss11 = A(:);  ss12 = B(:);  B1 = B(:,[1 4 7 2 5 8 3 6 9]);
ss21 = B1(:);
ii = [ii11; ii12; ii21];
jj = [jj11; jj12; jj21];
ss = [ss11; ss12; ss21];
kk = sparse(ii,jj,ss,2*N,2*N);

ff = accumarray(elem(:)+N, F(:),[2*N 1]);

% ------- Neumann boundary conditions -----------
elemN = bdStruct.elemD;
z1 = node(elemN(:,1),:); z2 = node(elemN(:,2),:);
e = z1-z2; % e = z2-z1
ne = [-e(:,2),e(:,1)]; % scaled ne
Du = pde.Du;
gradu1 = Du(z1); gradu2 = Du(z2);
F1 = sum(ne.*gradu1,2)./2; F2 = sum(ne.*gradu2,2)./2;
FN = [F1,F2];
ff = ff + accumarray(elemN(:), FN(:),[2*N 1]);

% --------- Dirichlet boundary conditions ---------------
eD = bdStruct.eD; g_D = pde.g_D;
id = eD+N;
isBdDof = false(2*N,1); isBdDof(id) = true;
bdDof = isBdDof; freeDof = (~isBdDof);
pD = node(eD,:);
U = zeros(2*N,1); U(bdDof) = g_D(pD);
ff = ff - kk*U;

% ------------ Solver -----------
U(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
u = U(N+1:end); w = U(1:N);
