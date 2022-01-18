function uh = varFEM1D(Th,pde,Vh,quadOrder)
%varFEM1D solves the two-point boundary value problem
%
%     -au'' + bu' + cu = f  in \Omega = (a,b), with
%     Dirichlet boundary condition u = g_D  on \Gamma_D = {a}, {b} or {a,b},
%     Neumann boundary condition   u = u' on \Gamma_N = {a,b} - \Gamma_D
%
% Copyright (C) Terence Yu.

% Quadrature orders for int1d and int2d
if nargin==2, Vh = 'P1'; quadOrder = 3; end % default: P1
if nargin==3, quadOrder = 3; end

%% Assemble stiffness matrix 
Coef  = {pde.a,   pde.b,   pde.c};
Test  = {'v.dx', 'v.val',  'v.val'};
Trial = {'u.dx',  'u.dx',  'u.val'};
kk = int1d(Th,Coef,Test,Trial,Vh,quadOrder);

%% Assemble the right hand side 
Coef = pde.f; Test = 'v.val';
ff = int1d(Th,Coef,Test,[],Vh,quadOrder);

%% Apply the Dirichlet boundary value conditions
uh = Applyboundary1D(Th,kk,ff,pde);