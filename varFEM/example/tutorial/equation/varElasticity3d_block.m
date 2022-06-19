function u = varElasticity3d_block(Th,pde,Vh,quadOrder)
%varElasticity3d_block Conforming Lagrange elements for linear elasticity
%equation in 3-D
% 
% Variational formulation based programming
%       u = [u1, u2, u3]
%       -div (sigma) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D, g3_D] on \Gamma_D
%       Neumann boundary condition   \sigma*n = g  on \Gamma_N
%       \sigma = (sigma_{ij}): stress tensor, 1<=i,j<=3
%

% Quadrature orders for int1d and int2d
if nargin==2, Vh = 'P1'; quadOrder = 4; end % default: P1
if nargin==3, quadOrder = 4; end

mu = pde.mu; lambda = pde.lambda; 

%% Get matrices of all pairs
% (vi.dx, uj.dx) 
Coef  = 1;  Test  = 'v.dx';  Trial = 'u.dx';  
Axx = assem3d(Th,Coef,Test,Trial,Vh,quadOrder);
% (vi.dx, uj.dy) 
Coef  = 1;  Test  = 'v.dx';  Trial = 'u.dy';  
Axy = assem3d(Th,Coef,Test,Trial,Vh,quadOrder);
% (vi.dx, uj.dz) 
Coef  = 1;  Test  = 'v.dx';  Trial = 'u.dz';  
Axz = assem3d(Th,Coef,Test,Trial,Vh,quadOrder);
% (vi.dy, uj.dx) 
Coef  = 1;  Test  = 'v.dy';  Trial = 'u.dx';  
Ayx = assem3d(Th,Coef,Test,Trial,Vh,quadOrder);
% (vi.dy, uj.dy) 
Coef  = 1;  Test  = 'v.dy'; Trial = 'u.dy';  
Ayy = assem3d(Th,Coef,Test,Trial,Vh,quadOrder);
% (vi.dy, uj.dz) 
Coef  = 1;  Test  = 'v.dy'; Trial = 'u.dz';  
Ayz = assem3d(Th,Coef,Test,Trial,Vh,quadOrder);
% (vi.dz, uj.dx) 
Coef  = 1;  Test  = 'v.dz'; Trial = 'u.dx';  
Azx = assem3d(Th,Coef,Test,Trial,Vh,quadOrder);
% (vi.dz, uj.dy) 
Coef  = 1;  Test  = 'v.dz'; Trial = 'u.dy';  
Azy = assem3d(Th,Coef,Test,Trial,Vh,quadOrder);
% (vi.dz, uj.dz) 
Coef  = 1;  Test  = 'v.dz'; Trial = 'u.dz';  
Azz = assem3d(Th,Coef,Test,Trial,Vh,quadOrder);

%% Get block stiffness matrix 
% (Eij(u):Eij(v))
A = [ Axx + 0.5*(Ayy+Azz)        0.5*Ayx           0.5*Azx;
           0.5*Axy         Ayy + 0.5*(Axx+Azz)     0.5*Azy;   
           0.5*Axz         0.5*Ayz           Azz + 0.5*(Axx+Ayy)   ];

% (div u,div v) 
B = [ Axx    Axy   Axz;
      Ayx    Ayy   Ayz
      Azx    Azy   Azz];

% stiffness matrix
kk = sparse(2*mu*A + lambda*B);

%% Assmeble right hand side
% F1
trf = eye(3); 
Coef = @(pz) pde.f(pz)*trf(:, 1);  Test = 'v.val';
F1 = assem3d(Th,Coef,Test,[],Vh,quadOrder);
% F2
Coef = @(pz) pde.f(pz)*trf(:, 2);  Test = 'v.val';
F2 = assem3d(Th,Coef,Test,[],Vh,quadOrder);
% F2
Coef = @(pz) pde.f(pz)*trf(:, 3);  Test = 'v.val';
F3 = assem3d(Th,Coef,Test,[],Vh,quadOrder);

%% Assemble Neumann boundary conditions
if ~isempty(Th.bdStr)
    Th.on = 1;
    g_N = pde.g_N; trg = eye(6);    
    g1 = @(p) g_N(p)*trg(:,[1,4,5]);   Cmat1 = interpFaceMat(g1,Th,quadOrder);
    g2 = @(p) g_N(p)*trg(:,[4,2,6]);   Cmat2 = interpFaceMat(g2,Th,quadOrder);
    g3 = @(p) g_N(p)*trg(:,[5,6,3]);   Cmat3 = interpFaceMat(g3,Th,quadOrder);

    Coef = Cmat1; Test = 'v.val';
    F1 = F1 + assem2d(Th,Coef,Test,[],Vh,quadOrder);
    Coef = Cmat2; Test = 'v.val';
    F2 = F2 + assem2d(Th,Coef,Test,[],Vh,quadOrder);
    Coef = Cmat3; Test = 'v.val';
    F3 = F3 + assem2d(Th,Coef,Test,[],Vh,quadOrder);
end

% load vector
ff = [F1; F2; F3];

%% Apply Dirichlet boundary conditions
on = 2 - 1*isempty(Th.bdStr);
g_D = pde.g_D;  
g_D1 = @(p) g_D(p)*[1;0;0];
g_D2 = @(p) g_D(p)*[0;1;0];
g_D3 = @(p) g_D(p)*[0;0;1];
gBc = {g_D1,g_D2,g_D3};
Vhvec = {Vh,Vh,Vh};
u = apply3d(on,Th,kk,ff,Vhvec,gBc); % note Vhvec