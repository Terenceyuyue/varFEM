function uh = FEM1D_variational(Th,pde,feSpace,quadOrder)

% Quadrature orders for int1d and int2d
if nargin==2, feSpace = 'P1'; quadOrder = 3; end % default: P1
if nargin==3, quadOrder = 3; end

% --------------- Stiffness matrix -----------------
Coef  = {pde.a,   pde.b,   pde.c};
Trial = {'u.dx',  'u.dx',  'u.val'};
Test  = {'v.dx', 'v.val',  'v.val'};
kk = int1d(Th,Coef,Trial,Test,feSpace,quadOrder);

% ------------------ Load vector -------------------
Coef = pde.f; Test = 'v.val';
ff = int1d(Th,Coef,[],Test,feSpace,quadOrder);

% --------- Boundary value conditions ---------
uh = Applyboundary1D(Th,kk,ff,pde);