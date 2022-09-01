function u = varElasticity_block(Th,pde,Vh,quadOrder)
%varElasticity_block Conforming P1 FEM of linear elasticity equation 
% variational formulation based programming in terms of block matrix
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

%% Get matrices of all pairs
% (vi.dx, uj.dx) 
cf = 1;
Coef  = cf;  Test  = 'v.dx';  Trial = 'u.dx';  
Axx = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
% (vi.dx, uj.dy) 
cf = 1;
Coef  = cf;  Test  = 'v.dx';  Trial = 'u.dy';  
Axy = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
% (vi.dy, uj.dx) 
cf = 1;
Coef  = cf;  Test  = 'v.dy';  Trial = 'u.dx';  
Ayx = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
% (vi.dy, uj.dy) 
cf = 1;
Coef  = cf;  Test  = 'v.dy'; Trial = 'u.dy';  
Ayy = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);

%% Get block stiffness matrix 
% (Eij(u):Eij(v))
A = [ Axx + 0.5*Ayy         0.5*Ayx;
           0.5*Axy         0.5*Axx + Ayy ];

% (div u,div v) 
B = [ Axx    Axy;
      Ayx    Ayy ];

% stiffness matrix
kk = sparse(2*mu*A + lambda*B);

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
if ~isempty(Th.bdStr)
    Th.on = 1;
    g_N = pde.g_N; trg = eye(3);
    
    g1 = @(p) g_N(p)*trg(:,[1,3]);
    Cmat1 = interpEdgeMat(g1,Th,quadOrder);
    Coef = Cmat1; Test = 'v.val';
    F1 = F1 + assem1d(Th,Coef,Test,[],Vh,quadOrder);
    
    g2 = @(p) g_N(p)*trg(:,[3,2]);
    Cmat2 = interpEdgeMat(g2,Th,quadOrder);
    Coef = Cmat2; Test = 'v.val';
    F2 = F2 + assem1d(Th,Coef,Test,[],Vh,quadOrder);
end

% load vector
ff = [F1; F2];

%% Apply Dirichlet boundary conditions
on = 2 - 1*isempty(Th.bdStr);
g_D = pde.g_D;  
g_D1 = @(p) g_D(p)*[1;0];
g_D2 = @(p) g_D(p)*[0;1];
gBc = {g_D1,g_D2};
Vhvec = {Vh,Vh};
u = apply2d(on,Th,kk,ff,Vhvec,gBc); % note Vhvec