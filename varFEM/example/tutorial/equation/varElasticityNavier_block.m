function u = varElasticityNavier_block(Th,pde,Vh,quadOrder)
%varElasticityNavier_block Conforming P1 elements discretization of linear elasticity equation of Navier form
% variational formulation based programming in terms of block matrix
%
%       u = [u1, u2]
%       -mu \Delta u - (lambda + mu)*grad(div(u)) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D.
%

% Quadrature orders for int1d and assem2d
if nargin==2, Vh = 'P1'; quadOrder = 3; end % default: P1
if nargin==3, quadOrder = 3; end

mu = pde.mu; lambda = pde.lambda; 

%% Assemble stiffness matrix
% (v1.grad, u1.grad), (v2.grad, u2.grad)
cf = 1;
Coef  = cf;  Test  = 'v.grad';  Trial = 'u.grad';  
A = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
% (v1.dx, u1.dx)
cf = 1;
Coef  = cf;  Test  = 'v.dx';  Trial = 'u.dx';  
B1 = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
% (v1.dx, u2.dy)
cf = 1;
Coef  = cf;  Test  = 'v.dx';  Trial = 'u.dy';  
B2 = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
% (v2.dy, u1.dx)
cf = 1;
Coef  = cf;  Test  = 'v.dy';  Trial = 'u.dx';  
B3 = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
% (v2.dy, u2.dy)
cf = 1;
Coef  = cf;  Test  = 'v.dy';  Trial = 'u.dy';  
B4 = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);

% kk
kk = [  mu*A+(lambda+mu)*B1,         (lambda+mu)*B2;
             (lambda+mu)*B3,    mu*A+(lambda+mu)*B4   ];
kk = sparse(kk);
    
%% Assemble right hand side
% F1
Coef = @(pz) pde.f(pz)*[1;0];  Test = 'v.val';
F1 = assem2d(Th,Coef,Test,[],Vh,quadOrder);

% F2
Coef = @(pz) pde.f(pz)*[0;1];  Test = 'v.val';
F2 = assem2d(Th,Coef,Test,[],Vh,quadOrder);

% F
ff = [F1; F2];

%% Apply Dirichlet boundary conditions
g_D = pde.g_D;  
on = 1;
g_D1 = @(p) g_D(p)*[1;0];
g_D2 = @(p) g_D(p)*[0;1];
gBc = {g_D1,g_D2};
Vhvec = {Vh,Vh};
u = apply2d(on,Th,kk,ff,Vhvec,gBc); % note Vhvec