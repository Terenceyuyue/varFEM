function u = elasticity2_variationalBlock(Th,pde,feSpace,quadOrder)
%Elasticity2_variationalBlock  Conforming P1 FEM of linear elasticity equation 
% variational formulation based programming in terms of block matrix
%       u = [u1, u2]
%       -div (sigma) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D
%       Neumann boundary condition   \sigma*n = g  on \Gamma_N
%       \sigma = (sigma_{ij}): stress tensor, 1<=i,j<=2

% Quadrature orders for int1d and int2d
if nargin==2, feSpace = 'P1'; quadOrder = 3; end % default: P1
if nargin==3, quadOrder = 3; end

% ------------------------ Mesh Th --------------------------
node = Th.node;  elem = Th.elem; bdStruct = Th.bdStruct;
Th.elem1D = bdStruct.elemN; Th.bdIndex1D = bdStruct.bdIndexN;
N = size(node,1); 
mu = pde.mu; lambda = pde.lambda; 

% ------------------ matrices of all pairs ------------------
% (vi.dx, uj.dx) 
cf = 1;
Coef  = {cf};  Trial = {'u.dx'};  Test  = {'v.dx'};
A1 = int2d(Th,Coef,Trial,Test,feSpace,quadOrder);
% (vi.dx, uj.dy) 
cf = 1;
Coef  = {cf};  Trial = {'u.dy'};  Test  = {'v.dx'};
A2 = int2d(Th,Coef,Trial,Test,feSpace,quadOrder);
% (vi.dy, uj.dx) 
cf = 1;
Coef  = {cf};  Trial = {'u.dx'};  Test  = {'v.dy'};
A3 = int2d(Th,Coef,Trial,Test,feSpace,quadOrder);
% (vi.dy, uj.dy) 
cf = 1;
Coef  = {cf};  Trial = {'u.dy'};  Test  = {'v.dy'};
A4 = int2d(Th,Coef,Trial,Test,feSpace,quadOrder);

% ---------------------- Stiffness matrix -----------------------
% (Eij(u):Eij(v))
A = [ A1 + 0.5*A4         0.5*A3;
           0.5*A2         0.5*A1 + A4 ];
A = 2*mu*A;

% (div u,div v) 
B = [ A1    A2;
      A3    A4 ];
B = lambda*B;

% stiffness matrix
kk = A + B;

% -------------------------- Load vector -------------------------
% F1
trf = eye(2); 
Coef = @(pz) pde.f(pz)*trf(:, 1);  Test = 'v.val';
F1 = int2d(Th,Coef,[],Test,feSpace,quadOrder);

% F2
trf = eye(2); 
Coef = @(pz) pde.f(pz)*trf(:, 2);  Test = 'v.val';
F2 = int2d(Th,Coef,[],Test,feSpace,quadOrder);

% ------------ Neumann boundary condition ----------------
if ~isempty(Th.elem1D)
    g_N = pde.g_N; trg = eye(3);
    
    g1 = @(p) g_N(p)*trg(:,[1,3]);
    Cmat1 = getMat1d(g1,Th,quadOrder);
    Coef = Cmat1; Test = 'v.val';
    F1 = F1 + int1d(Th,Coef,[],Test,feSpace,quadOrder);
    
    g2 = @(p) g_N(p)*trg(:,[3,2]);
    Cmat2 = getMat1d(g2,Th,quadOrder);
    Coef = Cmat2; Test = 'v.val';
    F2 = F2 + int1d(Th,Coef,[],Test,feSpace,quadOrder);
end

% load vector
ff = [F1; F2];

% ------------ Dirichlet boundary condition ----------------
g_D = pde.g_D;  eD = bdStruct.eD;
id = [eD; eD+N]; 
isBdNode = false(2*N,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
pD = node(eD,:);
u = zeros(2*N,1); uD = g_D(pD); u(bdDof) = uD(:);
ff = ff - kk*u;

% ------------------ Solver -------------------
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
