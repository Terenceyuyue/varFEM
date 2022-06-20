clc;clear;close all
%Poisson_VI_Uzawa solves the simplified friction problem using conforming
% finite element method in the lowest order case. The problem is iteratively 
% solved by using the Uzawa algorithm, hence the implementation is reduced to
% that of Poisson equation in each iteration.
%
% The problem is
%
%  D.E.
%     -\Delta u + cu = f,  in Omega,
%  Frictional boundary conditions
%     |grad(u)*n| <= g,  u*grad(u)*n + g*|u| = 0 on \Gamma_C
%  Dirichlet boundary conditions
%     u = g_D    on \Gamma_D
%  \partial(\Omega) = \Gamma_C + \Gamma_D
%
%   References
%   [1] F. Wang and H. Wei. Virtual element method for simplified friction problem. 
%       Appl. Math. Lett., 85: 125-131, 2018.
%   [2] B. Wu, F. Wang, and W. Han. Virtual element method for a frictional contact 
%       problem with normal compliance. Commun. Nonlinear Sci. Numer. Simul., 107: 
%       Paper No. 106125, 13 pp., 2022.
%

%% PDE data
c = 1e4; 
pde = PoissondataVI(c);
bdFriction = 'y==0'; % 1-string for frictional boundary condition

%% Mesh
[node,elem] = squaremesh([0 1 0 1],0.05);
% mesh info
Th = FeMesh2d(node,elem,bdFriction);
% Vh,quadOrder
Vh = 'P1'; quadOrder = 5;

%% Bilinear form
Coef = {1, c};
Test = {'v.grad', 'v.val'};
Trial = {'u.grad', 'u.val'};
kk = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);

%% Linear form
Coef = pde.f;
Test = 'v.val';
ff = assem2d(Th,Coef,Test,[],Vh,quadOrder);

%% Initialization of Uzawa iteration or Gamma_C is empty 
g_D = pde.g_D;  
on = 2 - 1*isempty(bdFriction); % 2-Dirichlet
uh0 = apply2d(on,Th,kk,ff,Vh,g_D);

%% Uzawa iteration
if ~isempty(bdFriction)
    Th.elem1d = Th.bdEdgeType{1}; % for friction
    Th.elem1dIdx = Th.bdEdgeIdxType{1};
    bdEdgeFri = Th.elem1d;
    g = pde.g_C;  
    z1 = node(bdEdgeFri(:,1),:); z2 = node(bdEdgeFri(:,2),:);
    gC = max([g(z1); g(z2)]);
    N = size(node,1);
    lambdah0 = ones(N,1); % Lagrange multiplizer   
     
    tol = 1e-8;  Err = 1;   
    iter = 0;  maxIt = 500;
    while Err>tol || iter<=maxIt      
        % assemble frictional boundary conditions
        Coef = gC*lambdah0;
        Test = 'v.val';
        ffNew = ff - assem1d(Th,Coef,Test,[],Vh,quadOrder);
        
        % apply Dirichlet boundary conditions
        uh = apply2d(on,Th,kk,ffNew,Vh,g_D);        
        
        % update the Lagrangian multiplier
        rho = 10;
        lambdah = lambdah0 + rho*gC*uh;
        lambdah = max(-ones(N,1), min(ones(N,1),lambdah));
        
        % compute errors between two steps
        Err = norm(uh-uh0);
        
        % update uh0 and lambdah0
        uh0 = uh; lambdah0 = lambdah;
        iter = iter + 1;
    end
end

%% The discrete solution
uh = uh0;
showresult(node,elem,pde.uexact,uh);

%% Interpolate to a 2D cartesian grid
figure, 
subplot(1,2,1),
varcontourf(node,elem,uh(1:Th.N),20);
title('Exact solution');
subplot(1,2,2),
ue = pde.uexact(node);
varcontourf(node,elem,ue(1:Th.N),20);
title('Numerical solution');