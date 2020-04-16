function [u,w] = biharmonicMixedFEM_variational(Th,pde,feSpace,quadOrder)
% BiharmonicMixedFEM_variational solves the biharmonic equation
% using the variational formulation based programming.
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

% The rate of convergence for u is optimal but for w is sub-optimal. For
% linear element, optimal order for w is also observed.

if nargin==2, feSpace = 'P1'; quadOrder = 3; end
if nargin==3, quadOrder = 3; end

% ------------------------ Mesh Th --------------------------
% node
node = Th.node;  N = size(node,1);

% -------------- Stiffness matrix -------------
Coef = { -1, 1, 1 }; 
Trial = {'u1.val', 'u2.grad', 'u1.grad'};
Test  = {'v1.val', 'v1.grad', 'v2.grad'};
kk = int2dvec(Th,Coef,Trial,Test,feSpace,quadOrder);

% -------------- Load vector -------------
Coef = pde.f; Test = 'v2.val';
ff = int2dvec(Th,Coef,[],Test,feSpace,quadOrder);

% ------- Neumann boundary conditions -----------
% 1D-mesh
bdStruct = Th.bdStruct;
Th.elem1D = bdStruct.elemD; Th.bdIndex1D = bdStruct.bdIndexD;
% Coef = @(p) pde.Du(p)*n;
Coef = getMat1d(pde.Du,Th,quadOrder);  Test = 'v1.val';
ff = ff + int1dvec(Th,Coef,[],Test,feSpace,quadOrder);

% --------- Dirichlet boundary conditions ---------------
stru = 'u2.val';
U = Applyboundary2Dvec(Th,kk,ff,pde,stru,feSpace);

U = reshape(U,[],2);
w = U(:,1);   u = U(:,2); 

