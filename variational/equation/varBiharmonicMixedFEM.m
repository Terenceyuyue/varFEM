function [u,w] = varBiharmonicMixedFEM(Th,pde,Vh,quadOrder)
%% VARBIHARMONICMIXEDFEM solves the biharmonic equation
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
% Get 1D mesh in boundary integrals 
elem = Th.elem; bdStruct = Th.bdStruct;
Th.elem1D = bdStruct.bdEdgeD; 
Th.bdType = 1;   % default = 2 for numerical integration on boundary ( see dof1d.m )
% Coef = @(p) pde.Du(p)*n;
Coef = getMat1d(pde.Du,Th,quadOrder);  Test = 'v1.val';
ff = ff + int1d(Th,Coef,Test,[],Vh,quadOrder);

%% Apply Dirichlet boundary conditions 
g_D = { [], pde.g_D };
U = Applyboundary2D(Th,kk,ff,g_D,Vh);
U = reshape(U,[],2);
w = U(:,1);   u = U(:,2); 

