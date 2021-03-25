function u = varElasticity_block(Th,pde,Vh,quadOrder)
%varElasticity_block Conforming P1 FEM of linear elasticity equation 
% variational formulation based programming in terms of block matrix
%       u = [u1, u2]
%       -div (sigma) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D
%       Neumann boundary condition   \sigma*n = g  on \Gamma_N
%       \sigma = (sigma_{ij}): stress tensor, 1<=i,j<=2
%
% Copyright (C) Terence Yu.

% Quadrature orders for int1d and int2d
if nargin==2, Vh = 'P1'; quadOrder = 3; end % default: P1
if nargin==3, quadOrder = 3; end

mu = pde.mu; lambda = pde.lambda; 

%% Get 1D mesh in boundary integrals
bdStruct = Th.bdStruct;
Th.elem1D = bdStruct.bdEdgeN; 

%% Get matrices of all pairs
% (vi.dx, uj.dx) 
cf = 1;
Coef  = {cf};  Test  = {'v.dx'};  Trial = {'u.dx'};  
A1 = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
% (vi.dx, uj.dy) 
cf = 1;
Coef  = {cf};  Test  = {'v.dx'};  Trial = {'u.dy'};  
A2 = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
% (vi.dy, uj.dx) 
cf = 1;
Coef  = {cf};  Test  = {'v.dy'};  Trial = {'u.dx'};  
A3 = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
% (vi.dy, uj.dy) 
cf = 1;
Coef  = {cf};  Test  = {'v.dy'}; Trial = {'u.dy'};  
A4 = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);

%% Get block stiffness matrix 
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

%% Assmeble right hand side
% F1
trf = eye(2); 
Coef = @(pz) pde.f(pz)*trf(:, 1);  Test = 'v.val';
F1 = assem2d(Th,Coef,Test,[],Vh,quadOrder);

% F2
trf = eye(2); 
Coef = @(pz) pde.f(pz)*trf(:, 2);  Test = 'v.val';
F2 = assem2d(Th,Coef,Test,[],Vh,quadOrder);

%% Assemble Neumann boundary conditions
if ~isempty(Th.elem1D)
    g_N = pde.g_N; trg = eye(3);
    
    g1 = @(p) g_N(p)*trg(:,[1,3]);
    Cmat1 = getMat1d(g1,Th,quadOrder);
    Coef = Cmat1; Test = 'v.val';
    F1 = F1 + assem1d(Th,Coef,Test,[],Vh,quadOrder);
    
    g2 = @(p) g_N(p)*trg(:,[3,2]);
    Cmat2 = getMat1d(g2,Th,quadOrder);
    Coef = Cmat2; Test = 'v.val';
    F2 = F2 + assem1d(Th,Coef,Test,[],Vh,quadOrder);
end

% load vector
ff = [F1; F2];

%% Apply Dirichlet boundary conditions
fixedNode = bdStruct.bdNodeIdx;
g_D = pde.g_D;  N = Th.N;
id = [fixedNode; fixedNode+N]; 
isBdNode = false(2*N,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
nodeD = Th.node(fixedNode,:);
u = zeros(2*N,1); uD = g_D(nodeD); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Direct Solver
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
