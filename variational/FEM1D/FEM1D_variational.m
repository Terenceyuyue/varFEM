function uh = FEM1D_variational(Th,pde,Vh,quadOrder)

% Quadrature orders for int1d and int2d
if nargin==2, Vh = 'P1'; quadOrder = 3; end % default: P1
if nargin==3, quadOrder = 3; end

% --------------- Stiffness matrix -----------------
Coef  = {pde.a,   pde.b,   pde.c};
Test  = {'v.dx', 'v.val',  'v.val'};
Trial = {'u.dx',  'u.dx',  'u.val'};
kk = int1d(Th,Coef,Test,Trial,Vh,quadOrder);

% ------------------ Load vector -------------------
Coef = pde.f; Test = 'v.val';
ff = int1d(Th,Coef,Test,[],Vh,quadOrder);

% --------- Boundary value conditions ---------
uh = Applyboundary1D(Th,kk,ff,pde);