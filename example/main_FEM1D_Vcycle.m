clc;clear;close all
% ------------ Initial mesh and boundary conditions ----------
a = 0; b = 1;
N0 = 3;  % number of initial nodes
x = linspace(a,b,N0)'; % uniform
%x = [a; sort(rand(N0-2,1)); b]; % random
node = x; elem = [(1:N0-1)', (2:N0)']; % initial mesh
J = 3; % solve it directly when J<=1

Neumann = [1]; Dirichlet = [N0];
bdStruct = struct('Dirichlet', Dirichlet, 'Neumann', Neumann);

% ------------- Final mesh and prolongation matrices -----------
[node,elem,Pro,Res] = Mesh1DVcycle(node,elem,J);

% ------------------- PDE -----------------
acoef = -1;  bcoef = 0;  ccoef = 0;
para = struct('acoef',acoef, 'bcoef',bcoef, 'ccoef',ccoef);
pde = pdedata1D(para);

% ----------- FEM1DVcycle -----------
u = FEM1DVcycle(node,elem,pde,bdStruct,Pro,Res);

% -------- error analysis -----------
[node,id] = sort(node);
u = u(id);
uexact = pde.uexact(node);
figure,plot(node,u,'k-',node,uexact,'r-*','linewidth',1);
xlabel('x'); ylabel('u');
legend('Numerical solution','Exact solution')
Err = u-uexact;
figure, plot(node,Err,'k','linewidth',1); legend('Absolute error');

