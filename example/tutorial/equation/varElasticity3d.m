function u = varElasticity3d(Th,pde,Vh,quadOrder)
%varElasticity3d Conforming Lagrange elements for linear elasticity
%equation in 3-D
% 
% Variational formulation based programming
%       u = [u1, u2, u3]
%       -div (sigma) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D, g3_D] on \Gamma_D
%       Neumann boundary condition   \sigma*n = g  on \Gamma_N
%       \sigma = (sigma_{ij}): stress tensor, 1<=i,j<=3
%

mu = pde.mu; lambda = pde.lambda; 

%% Assemble stiffness matrix 
% (Eij(u):Eij(v))
Coef = { 1, 1, 1, 0.5, 0.5, 0.5 }; 
Test  = {'v1.dx', 'v2.dy', 'v3.dz', ...
    'v1.dy + v2.dx', 'v1.dz + v3.dx', 'v2.dz + v3.dy'};
Trial = {'u1.dx', 'u2.dy', 'u3.dz', ...
    'u1.dy + u2.dx', 'u1.dz + u3.dx', 'u2.dz + u3.dy'};
A = int3d(Th,Coef,Test,Trial,Vh,quadOrder);

% (div u,div v) 
Coef = 1; 
Test  = 'v1.dx + v2.dy + v3.dz' ;
Trial = 'u1.dx + u2.dy + u3.dz';
B = int3d(Th,Coef,Test,Trial,Vh,quadOrder);

% stiffness matrix
kk = 2*mu*A + lambda*B;

%% Assemble the right hand side 
Coef = pde.f;  Test = 'v.val';
ff = int3d(Th,Coef,Test,[],Vh,quadOrder);

%% Assemble Neumann boundary conditions
if ~isempty(Th.bdStr)
    Th.on = 1;  

    g_N = pde.g_N; trg = eye(6);    
    g1 = @(p) g_N(p)*trg(:,[1,4,5]);   Cmat1 = interpFaceMat(g1,Th,quadOrder);
    g2 = @(p) g_N(p)*trg(:,[4,2,6]);   Cmat2 = interpFaceMat(g2,Th,quadOrder);
    g3 = @(p) g_N(p)*trg(:,[5,6,3]);   Cmat3 = interpFaceMat(g3,Th,quadOrder);
    
    Coef = {Cmat1, Cmat2, Cmat3};  Test = 'v.val';  
    ff = ff + int2d(Th,Coef,Test,[],Vh,quadOrder);
end

%% Apply Dirichlet boundary conditions
on = 2 - 1*isempty(Th.bdStr);
g_D1 = @(p) pde.g_D(p)*[1;0;0];
g_D2 = @(p) pde.g_D(p)*[0;1;0];
g_D3 = @(p) pde.g_D(p)*[0;0;1];
g_D = {g_D1, g_D2, g_D3};
u = apply3d(on,Th,kk,ff,Vh,g_D);