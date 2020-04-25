function [u,w] = biharmonicMixedFEM_variationalBlock(Th,pde)
% BiharmonicMixedFEM_variationalBlock solves the biharmonic equation using
% the mixed FEM
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

node = Th.node; bdStruct = Th.bdStruct;
N = size(node,1);

% -------------- Stiffness matrix -------------
% matrix A 
cf = @(p) 1 + 0*p(:,1); 
Coef = {cf}; Test = {'v.val'}; Trial = {'u.val'}; 
A = -assem2d(Th,Coef,Test,Trial);
% matrix B
cf = @(p) 1 + 0*p(:,1); 
Coef = {cf}; Test = {'v.grad'}; Trial = {'u.grad'}; 
B = assem2d(Th,Coef,Test,Trial);
% kk 
O = zeros(size(B));
kk = [A,  B;  B', O];

% -------------- Load vector -------------
Coef = pde.f; Test = 'v.val';
ff = assem2d(Th,Coef,Test,[]);
O = zeros(size(ff));
ff = [O; ff];

% ------- Neumann boundary conditions -----------
Th.elem1D = bdStruct.elemD;
%Coef = @(p) pde.Du(p)*n;
Coef = getMat1d(pde.Du,Th);
ff(1:N) = ff(1:N) + assem1d(Th,Coef,Test,[]);

% --------- Dirichlet boundary conditions ---------------
eD = bdStruct.eD; g_D = pde.g_D;
id = eD+N;
isBdDof = false(2*N,1); isBdDof(id) = true;
bdDof = isBdDof; freeDof = (~isBdDof);
pD = node(eD,:);
U = zeros(2*N,1); U(bdDof) = g_D(pD);
ff = ff - kk*U;

% ------------ Solver -----------
U(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
w = U(1:N);  u = U(N+1:end); 

