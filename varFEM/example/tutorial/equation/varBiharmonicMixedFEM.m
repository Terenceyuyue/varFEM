function [u,w] = varBiharmonicMixedFEM(Th,pde,Vh,quadOrder)
%varBiharmonicMixedFEM solves the biharmonic equation using the variational 
% formulation based programming
%     Laplace^2 u = f;   
%               u = g_D;
%            Dn u = g_N.
%
% by writing in a mixed form
%
%       - Laplace u = w
%       - Laplace w = f
%
% Unlike the conforming or nonconforming FEMs, the second boundary condition
% is imposed as a Neumann boundary condition.

% The rate of convergence for u is optimal but for w is sub-optimal:
% - For linear element, optimal order for w is also observed.
% - For P2 element, L2(w) is 1.5, H1(w) is 0.5
% - For P3 element, L2(w) is 2.5, H1(w) is 1.5
% For P2 and P3, the rate of w has the behaviour of Laplace(u)
%

if nargin==2, Vh = repmat({'P1'},1,2); quadOrder = 3; end
if nargin==3, quadOrder = 3; end

%% Assemble stiffness matrix
Coef = { -1, 1, 1 }; 
Test  = {'v1.val', 'v1.grad', 'v2.grad'};
Trial = {'u1.val', 'u2.grad', 'u1.grad'};
kk = int2d(Th,Coef,Test,Trial,Vh,quadOrder);

%% Assemble right hand side
Coef = pde.f; Test = 'v2.val';
ff = int2d(Th,Coef,Test,[],Vh,quadOrder);

%% Assemble Neumann boundary conditions
% Get 1D mesh for boundary integrals 
Th.on = 1;
% Th.elem1d = Th.bdEdgeType{1}; 
% Th.elem1dIdx = Th.bdEdgeIdxType{1};

% Coef = @(p) pde.Du(p)*n;
Coef = interpEdgeMat(pde.Du,Th,quadOrder);  
Test = 'v1.val';
ff = ff + int1d(Th,Coef,Test,[],Vh,quadOrder);

%% Apply Dirichlet boundary conditions 
g_D = { [], pde.g_D };
on = 1;
U = apply2d(on,Th,kk,ff,Vh,g_D);
U = reshape(U,[],2);
w = U(:,1);   u = U(:,2); 