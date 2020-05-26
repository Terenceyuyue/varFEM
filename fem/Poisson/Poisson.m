function u = Poisson(node,elem,pde,bdStruct)
% Poisson solves Poisson equation with P1 linear element (2D).
%
%    - Delta u = f   in \Omega, with
%    Dirichlet boundary conditions u = g_D on \Gamma_D,
%    Neumann boundary conditions   grad(u)*n = g_N on \Gamma_N.

N = size(node,1); NT = size(elem,1); Ndof = 3;
f = pde.f;
% -------------- Sparse assembling index ----------------
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

% -------------- Assemble stiffness matrix -------------
[Dphi,area] = gradbasis(node,elem);
K = zeros(NT,Ndof^2); % straighten
s = 1;
for i = 1:Ndof
    for j = 1:Ndof
        K(:,s) = sum(Dphi(:,:,i).*Dphi(:,:,j),2).*area;
        s = s+1;
    end
end
kk = sparse(ii,jj,K(:),N,N);

% ------------- Assemble load vector ------------
% Gauss quadrature rule
[lambda,weight] = quadpts(2);
F = zeros(NT,3); % straighten
for p = 1:length(weight)
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    fxy = f(pxy);  fv = fxy*lambda(p,:); % [f*phi1, f*phi2, f*phi3] at (xp,yp)
    F = F + weight(p)*fv;
end
F = repmat(area,1,3).*F;  % F = area.*F;
ff = accumarray(elem(:), F(:),[N 1]);

% ------- Neumann boundary conditions -----------
elemN = bdStruct.elemN;
if ~isempty(elemN)
    z1 = node(elemN(:,1),:); z2 = node(elemN(:,2),:);
    e = z1-z2; % e = z2-z1
    ne = [-e(:,2),e(:,1)]; % scaled ne
    Du = pde.Du;
    gradu1 = Du(z1); gradu2 = Du(z2);
    F1 = sum(ne.*gradu1,2)./2; F2 = sum(ne.*gradu2,2)./2;
    FN = [F1,F2];
    ff = ff + accumarray(elemN(:), FN(:),[N 1]);
end

% --------- Dirichlet boundary conditions ---------------
eD = bdStruct.eD; g_D = pde.g_D;
isBdNode = false(N,1); isBdNode(eD) = true;
bdNode = (isBdNode); freeNode = (~isBdNode);
pD = node(bdNode,:);
u = zeros(N,1); u(bdNode) = g_D(pD);
ff = ff - kk*u;

% ------------ Solver -----------
u(freeNode) = kk(freeNode,freeNode)\ff(freeNode);

