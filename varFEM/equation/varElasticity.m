function u = varElasticity(Th,pde,Vh,quadOrder)
%varElasticity Conforming Lagrange elements for linear elasticity equation 
% 
% Variational formulation based programming
%       u = [u1, u2]
%       -div (sigma) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D
%       Neumann boundary condition   \sigma*n = g  on \Gamma_N
%       \sigma = (sigma_{ij}): stress tensor, 1<=i,j<=2
%

% Quadrature orders for int1d and int2d
if nargin==2, Vh = 'P1'; quadOrder = 3; end % default: P1
if nargin==3, quadOrder = 3; end

mu = pde.mu; lambda = pde.lambda; 

%% Get 1D mesh for boundary integrals
bdStr = Th.bdStr;
if ~isempty(bdStr)
    Th.elem1d = Th.bdEdgeType{1};
    Th.bdEdgeIdx1d = Th.bdEdgeIdxType{1};
end

%% Assemble stiffness matrix 
% (Eij(u):Eij(v))
Coef = { 1, 1, 0.5 }; 
Test  = {'v1.dx', 'v2.dy', 'v1.dy + v2.dx'};
Trial = {'u1.dx', 'u2.dy', 'u1.dy + u2.dx'};
A = int2d(Th,Coef,Test,Trial,Vh,quadOrder);
A = 2*mu*A;

% (div u,div v) 
Coef = { 1 }; 
Test  = { 'v1.dx + v2.dy' };
Trial = { 'u1.dx + u2.dy' };
B = int2d(Th,Coef,Test,Trial,Vh,quadOrder);
B = lambda*B;

% stiffness matrix
kk = A + B;

%% Assemble the right hand side 
Coef = pde.f;  Test = 'v.val';
ff = int2d(Th,Coef,Test,[],Vh,quadOrder);

%% Assemble Neumann boundary conditions
if ~isempty(bdStr)
    g_N = pde.g_N; trg = eye(3);
    
    g1 = @(p) g_N(p)*trg(:,[1,3]);   Cmat1 = getMatEdge(g1,Th,quadOrder);
    g2 = @(p) g_N(p)*trg(:,[3,2]);   Cmat2 = getMatEdge(g2,Th,quadOrder);
        
    Coef = {Cmat1, Cmat2};  Test = 'v.val';    
    ff = ff + int1d(Th,Coef,Test,[],Vh,quadOrder);
end

%% Apply Dirichlet boundary conditions
on = 2 - 1*isempty(bdStr);
tru = eye(2);
g_D1 = @(pz) pde.g_D(pz)*tru(:, 1);
g_D2 = @(pz) pde.g_D(pz)*tru(:, 2);
g_D = {g_D1, g_D2};
u = apply2d(on,Th,kk,ff,g_D,Vh);