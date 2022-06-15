function [u,w] = varBiharmonicMixedFEM_block(Th,pde)
%varBiharmonicMixedFEM_block solves the biharmonic equation using the mixed FEM
% variational formulation based programming in terms of block matrix
%
%     Laplace^2 u = f;   [a1,b1] * [a2,b2]
%               u = g_D;
%            Dn u = g_N.
%
% by writing in a mixed form
%
%       - Laplace u = w
%       - Laplace w = f
%
% Unlike the conforming or nonconforming FEMs, the second boundary condition
% is an Neumann boundary condition in this case.
%

node = Th.node; N = size(node,1);
quadOrder = 5;

%% Assemble stiffness matrix
% matrix A 
Coef = 1;  Test = 'v.val';  Trial = 'u.val'; 
A = -assem2d(Th,Coef,Test,Trial,'P1',quadOrder);
% matrix B
Coef = 1;  Test = 'v.grad';  Trial = 'u.grad'; 
B = assem2d(Th,Coef,Test,Trial,'P1',quadOrder);
% kk 
O = zeros(size(B));
kk = [A,  B;  B', O];
kk = sparse(kk);

%% Assemble right-hand side
Coef = pde.f;  Test = 'v.val';
ff = assem2d(Th,Coef,Test,[],'P1',quadOrder);
O = zeros(size(ff));
ff = [O; ff];

%% Assemble Neumann boundary conditions
Th.elem1d = Th.bdEdge; % all boundary edges
%Th.bdEdgeIdx1 = Th.bdEdgeIdx;
%Coef = @(p) pde.Du(p)*n;
Coef = interpEdgeMat(pde.Du,Th,quadOrder);
Test = 'v.val';
ff(1:N) = ff(1:N) + assem1d(Th,Coef,Test,[],'P1',quadOrder);

%% Apply Dirichlet boundary conditions
on = 1;
g_D = pde.g_D;
gBc = {[],g_D};
Vhvec = {'P1','P1'};
U = apply2d(on,Th,kk,ff,Vhvec,gBc); % note Vhvec
w = U(1:N);  u = U(N+1:end); 