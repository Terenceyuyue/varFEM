function u = elasticity3_variationalBlock(Th,pde,Vh,quadOrder)
%Elasticity3_variationalBlock  Conforming P1 elements discretization of linear elasticity equation
% variational formulation based programming in terms of block matrix
%
%       u = [u1, u2]
%       -mu \Delta u - (lambda + mu)*grad(div(u)) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D.
%


% Quadrature orders for int1d and assem2d
if nargin==2, Vh = 'P1'; quadOrder = 3; end % default: P1
if nargin==3, quadOrder = 3; end

% ------------------------ Mesh Th --------------------------
node = Th.node; bdStruct = Th.bdStruct; N = size(node,1); 
mu = pde.mu; lambda = pde.lambda; 

% ---------------------- Stiffness matrix -----------------------
% (v1.grad, u1.grad), (v2.grad, u2.grad)
cf = 1;
Coef  = {cf};  Test  = {'v.grad'};  Trial = {'u.grad'};  
A = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);

% (v1.dx, u1.dx) and (v2.dx, u1.dx)
cf = 1;
Coef  = {cf};  Test  = {'v.dx'};  Trial = {'u.dx'};  
B1 = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);

% (v1.dx, u2.dy) and (v2.dx, u2.dy)
cf = 1;
Coef  = {cf};  Test  = {'v.dx'};  Trial = {'u.dy'};  
B2 = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);

% kk
kk = [  mu*A+(lambda+mu)*B1,         (lambda+mu)*B2;
             (lambda+mu)*B1,    mu*A+(lambda+mu)*B2   ];
    
% -------------------------- Load vector -------------------------
% F1
trf = eye(2); 
Coef = @(pz) pde.f(pz)*trf(:, 1);  Test = 'v.val';
F1 = assem2d(Th,Coef,Test,[],Vh,quadOrder);

% F2
trf = eye(2); 
Coef = @(pz) pde.f(pz)*trf(:, 2);  Test = 'v.val';
F2 = assem2d(Th,Coef,Test,[],Vh,quadOrder);

% F
ff = [F1; F2];

% ------------ Dirichlet boundary condition ----------------
g_D = pde.g_D;  eD = bdStruct.eD;
id = [eD; eD+N]; 
isBdNode = false(2*N,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
pD = node(eD,:);
u = zeros(2*N,1); uD = g_D(pD); u(bdDof) = uD(:);
ff = ff - kk*u;

% ------------------ Solver -------------------
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
