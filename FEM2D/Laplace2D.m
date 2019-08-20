% Laplace2D: solve Poisson equation with P1 linear element.
%
%    - Delta u = f   in \Omega, with
%    Dirichlet boundary conditions u = g_D on \Gamma_D,
%    Neumann boundary conditions   grad(u)*n = g_N on \Gamma_N.

clc;clear;close all;
%  -------------- Mesh and boundary conditions --------------
[node,elem,eD,elemN] = meshLap();
N = size(node,1);

% ------------------------- PDE data ------------------------
% pde = struct('uexact',@uexact, 'f',@f, 'g_N',@g_N,  'g_D',@g_D);
pde = pdedata();

% -------------- Stiff matrix and load vector -------------
[kk,ff] = Get_kk_ff(node,elem,pde);

% ------- Neumann boundary conditions -----------
z1 = node(elemN(:,1),:); % starting points
z2 = node(elemN(:,2),:); % ending points
h = sqrt(sum((z2-z1).^2,2));
g_N = pde.g_N;
F1 = h./2.*g_N(z1); F2 = h./2.*g_N(z2);
FN = [F1,F2];
ff = ff + accumarray(elemN(:), FN(:),[N 1]);


% --------- Dirichlet boundary conditions ---------------
g_D = pde.g_D;
N = size(node,1);
isBdNode = false(N,1); isBdNode(eD) = true;
bdNode = find(isBdNode); freeNode = find(~isBdNode);
pD = node(bdNode,:);
u = zeros(N,1); u(bdNode) = g_D(pD);
ff = ff - kk*u;


% ------------ Solve the linear system -----------
u(freeNode) = kk(freeNode,freeNode)\ff(freeNode);

% --------- error analysis ---------------
uexact = pde.uexact;
ue = uexact(node);
format shorte
Err = norm(u-ue)/norm(ue)

figure,showsolution(node,elem,u);
figure,showsolution(node,elem,ue);

Eabs = u-ue;  % Absolute errors
figure,showsolution(node,elem,Eabs);
zlim('auto');

