function u = varElasticityNavierCR_block(Th,pde,Vh,quadOrder)
%varElasticityNavierCR_block CR element discretization of linear elasticity equation of Navier form
% variational formulation based programming in terms of block matrix
%
%       u = [u1, u2]
%       -mu \Delta u - (lambda + mu)*grad(div(u)) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D.
%
% Copyright (C) Terence Yu.

% Quadrature orders for int1d and assem2d
if nargin==2, Vh = 'CR'; quadOrder = 3; end 
if nargin==3, quadOrder = 3; end

mu = pde.mu; lambda = pde.lambda; 

%% Assemble stiffness matrix
% (v1.grad, u1.grad), (v2.grad, u2.grad)
cf = 1;
Coef  = {cf};  Test  = {'v.grad'};  Trial = {'u.grad'};  
A = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
% (v1.dx, u1.dx)
cf = 1;
Coef  = {cf};  Test  = {'v.dx'};  Trial = {'u.dx'};  
B1 = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
% (v1.dx, u2.dy)
cf = 1;
Coef  = {cf};  Test  = {'v.dx'};  Trial = {'u.dy'};  
B2 = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
% (v2.dy, u1.dx)
cf = 1;
Coef  = {cf};  Test  = {'v.dy'};  Trial = {'u.dx'};  
B3 = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
% (v2.dy, u2.dy)
cf = 1;
Coef  = {cf};  Test  = {'v.dy'};  Trial = {'u.dy'};  
B4 = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);

% kk
kk = [  mu*A+(lambda+mu)*B1,         (lambda+mu)*B2;
             (lambda+mu)*B3,    mu*A+(lambda+mu)*B4   ];
kk = sparse(kk);
    
%% Assemble right hand side
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

%% Apply Dirichlet boundary conditions
node = Th.node; bdStruct = Th.bdStruct;
NE = Th.NE;
isBdDof = false(2*NE,1);
bdEdgeIdxD = bdStruct.bdEdgeIdxD;
fixedNode = [bdEdgeIdxD; bdEdgeIdxD+NE];
isBdDof(fixedNode) = true;
bdDof = (isBdDof); freeDof = (~isBdDof);
bdEdgeD = bdStruct.bdEdgeD;
zm = (node(bdEdgeD(:,1),:)+node(bdEdgeD(:,2),:))/2; 
bdval = pde.g_D(zm);
u = zeros(2*NE,1); u(bdDof) = bdval(:);
ff = ff - kk*u;

%% Direct solver 
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
