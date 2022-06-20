function [uh,ph] = varStokes_block(Th,pde,Vh,quadOrder)
%varStokes_block Taylor-Hood (P2-P1) FEM for Stokes problem
% variational formulation based programming in terms of block matrix
%       u = [u1, u2]
%       -div(mu*grad u) + grad p = f in \Omega,
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma
%

% Quadrature orders for int1d and int2d
if nargin==2, Vh = {'P2','P2','P1'}; quadOrder = 5; end % default: Taylor-Hood
if nargin==3, quadOrder = 4; end

%% Assemble stiffness matrix
% A1,A2
Coef = 1;  Test = 'v1.grad';  Trial = 'u1.grad';  
A1 = assem2d(Th,Coef,Test,Trial,Vh(1),quadOrder);
A2 = A1;

% B1,B2
Coef = -1;  Test = 'v1.dx'; Trial = 'p.val';  
B1 = assem2d(Th,Coef,Test,Trial,Vh([1,3]),quadOrder);
Coef = -1;  Test = 'v2.dy'; Trial = 'p.val';  
B2 = assem2d(Th,Coef,Test,Trial, Vh([1,3]), quadOrder); 

% C1,C2
C1 = B1'; C2 = B2';

% D
eps = 1e-10;
Coef = -eps;  Test = 'p.val'; Trial = 'q.val';  
D = assem2d(Th,Coef,Test,Trial,Vh(3),quadOrder);

% stiffness matrix
NNdof = size(A1,1);
O = zeros(NNdof,NNdof); 
kk = [ A1,   O,   B1;
       O,    A2,  B2;
       C1,   C2,  D  ];
kk = sparse(kk);

%% Assemble right hand side 
ff = zeros(size(kk,1),1);
% F1
trf = eye(2); 
Coef = @(pz) pde.f(pz)*trf(:, 1);  Test = 'v1.val';
F1 = assem2d(Th,Coef,Test,[],Vh(1),quadOrder);

% F2
trf = eye(2); 
Coef = @(pz) pde.f(pz)*trf(:, 2);  Test = 'v2.val';
F2 = assem2d(Th,Coef,Test,[],Vh(2),quadOrder);

ff(1:2*NNdof) = [F1;F2];

%% Apply Dirichlet boundary conditions
tru = eye(2);
g_D1 = @(pz) pde.g_D(pz)*tru(:, 1);
g_D2 = @(pz) pde.g_D(pz)*tru(:, 2);
g_D = {g_D1, g_D2, []};
on = 1;
U = apply2d(on,Th,kk,ff,Vh,g_D);
uh = [U(1:NNdof), U(NNdof+1:2*NNdof)]; % uh = [u1h, u2h]
ph = U(2*NNdof+1:end);