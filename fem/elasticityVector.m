function u = elasticityVector(node,elem,pde,bdStruct)
%ElasticityVector  Conforming P1 FEM of linear elasticity equation
% Programming in the vectorized finite element space
%
%       u = [u1, u2]
%       -div (sigma) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D
%       Neumann boundary condition   \sigma*n = g  on \Gamma_N
%       \sigma = (sigma_{ij}): stress tensor, 1<=i,j<=2

N = size(node,1); NT = size(elem,1); Ndof = 3; Ndof2 = 2*Ndof;
mu = pde.mu; lambda = pde.lambda; f = pde.f;
% -------------- Compute (Dibase,Djbase) --------------------
[Dphi,area] = gradbasis(node,elem);
Dbase = cell(3,3);
for i = 1:2
    for j = 1:2
        Dbase{1,1}{i,j} = Dphi(:,i,1).*Dphi(:,j,1).*area; % b_st ---> (phi_s, phi_t)
        Dbase{1,2}{i,j} = Dphi(:,i,1).*Dphi(:,j,2).*area;
        Dbase{1,3}{i,j} = Dphi(:,i,1).*Dphi(:,j,3).*area;
        Dbase{2,1}{i,j} = Dphi(:,i,2).*Dphi(:,j,1).*area;
        Dbase{2,2}{i,j} = Dphi(:,i,2).*Dphi(:,j,2).*area;
        Dbase{2,3}{i,j} = Dphi(:,i,2).*Dphi(:,j,3).*area;
        Dbase{3,1}{i,j} = Dphi(:,i,3).*Dphi(:,j,1).*area;
        Dbase{3,2}{i,j} = Dphi(:,i,3).*Dphi(:,j,2).*area;
        Dbase{3,3}{i,j} = Dphi(:,i,3).*Dphi(:,j,3).*area;
    end
end

% -------- Sparse assembling indices -----------
nnz = NT*Ndof2^2;
ii = zeros(nnz,1); jj = zeros(nnz,1);
id = 0;
for i = 1:Ndof2
    for j = 1:Ndof2
        ix = (i==1 || i==3 || i==5);  jx = (j==1 || j==3 || j==5);
        i1 = (i+ix)/2; j1 = (j+jx)/2;
        ii(id+1:id+NT) = 2*elem(:,i1)-ix; % zi        
        jj(id+1:id+NT) = 2*elem(:,j1)-jx; % zj            
        id = id + NT;
    end
end

% ----------- Assemble stiffness matrix -----------
K = zeros(NT,Ndof2^2);
s1 = 1; s2 = s1+6; 
for s = 1:3
    for t = 1:3
        k11 = (lambda+2*mu)*Dbase{s,t}{1,1} + mu*Dbase{s,t}{2,2};
        k12 = lambda*Dbase{s,t}{1,2} + mu*Dbase{s,t}{2,1};
        k21 = lambda*Dbase{s,t}{2,1} + mu*Dbase{s,t}{1,2};
        k22 = (lambda+2*mu)*Dbase{s,t}{2,2} + mu*Dbase{s,t}{1,1};    
        K(:,s1) = k11; K(:,s1+1) = k12; s1 = s1+2;
        K(:,s2) = k21; K(:,s2+1) = k22; s2 = s2+2;
    end
    s1 = s1+Ndof2; s2 = s2+Ndof2;
end
kk = sparse(ii,jj,K(:),2*N,2*N);

% ------------- Assemble load vector ------------
% % mid-point quadrature rule
% x1 = node(elem(:,1),1); y1 = node(elem(:,1),2);
% x2 = node(elem(:,2),1); y2 = node(elem(:,2),2);
% x3 = node(elem(:,3),1); y3 = node(elem(:,3),2);
% xc = 1/3*(x1+x2+x3); yc = 1/3*(y1+y2+y3); pc = [xc,yc];
% F1 = f(pc).*area./3; F2 = F1; F3 = F1;
% F = [F1,F2,F3];

% Gauss quadrature rule
[lambda,weight] = quadpts(2);
F = zeros(NT,Ndof2); 
weight = [weight(:),weight(:)];
for iel = 1:NT
    vK = node(elem(iel,:),:); % vertices of K
    xy = lambda*vK;  fxy = f(xy); % fxy = [f1xy,f2xy]
    fv1 = fxy.*[lambda(:,1),lambda(:,1)]; % (f,phi1)
    fv2 = fxy.*[lambda(:,2),lambda(:,2)]; % (f,phi2)
    fv3 = fxy.*[lambda(:,3),lambda(:,3)]; % (f,phi3)
    
    F(iel,1:2) = area(iel)*dot(weight,fv1);
    F(iel,3:4) = area(iel)*dot(weight,fv2);
    F(iel,5:6) = area(iel)*dot(weight,fv3);
end
subs = zeros(NT,Ndof2);
subs(:,[1,3,5]) = 2*elem-1; subs(:,[2,4,6]) = 2*elem;
ff = accumarray(subs(:),F(:),[2*N,1]);

% ------------ Neumann boundary condition ----------------
elemN = bdStruct.elemN;
if ~isempty(elemN)
    g_N = pde.g_N;
    z1 = node(elemN(:,1),:); z2 = node(elemN(:,2),:);
    e = z1-z2;  % e = z2-z1
    ne = [-e(:,2),e(:,1)]; % scaled ne
    Sig1 = g_N(z1); Sig2 = g_N(z2);
    F11 = sum(ne.*Sig1(:,[1,3]),2)./2; F12 = sum(ne.*Sig2(:,[1,3]),2)./2;
    F21 = sum(ne.*Sig1(:,[3,2]),2)./2; F22 = sum(ne.*Sig2(:,[3,2]),2)./2;
    FN = [F11,F12,F21,F22];
    subs = [2*elemN(:)-1; 2*elemN(:)];
    ff = ff + accumarray(subs, FN(:),[2*N 1]);
end

% ------------ Dirichlet boundary condition ----------------
g_D = pde.g_D;  eD = bdStruct.eD;
isBdNode = false(N,1); isBdNode(eD) = true;
bdNode = find(isBdNode); freeNode = find(~isBdNode);
pD = node(bdNode,:);
bdDof = [2*bdNode-1; 2*bdNode]; freeDof = [2*freeNode-1;2*freeNode];
u = zeros(2*N,1); uD = g_D(pD); u(bdDof) = uD(:);
ff = ff - kk*u;

% ------------------ Solver -------------------
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);