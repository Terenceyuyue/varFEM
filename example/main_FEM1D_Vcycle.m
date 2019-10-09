clc;clear;close all
% ------------ Initial mesh and boundary conditions ----------
a = 0; b = 1;
N0 = 3;  x = linspace(a,b,N0)'; % N0: number of initial nodes
node0 = x; elem0 = [(1:N0-1)', (2:N0)']; % initial mesh
J = 4; % level length, J>=2

Neumann = 1; Dirichlet = N0;
bdStruct = struct('Dirichlet', Dirichlet, 'Neumann', Neumann);

% ------------- mesh and prolongation -----------
[node,elem,Pro,Res] = MeshVcycle(node0,elem0,J);

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
figure,plot(node,u,'k',node,uexact,'r*');
xlabel('x'); ylabel('u');
legend('Numerical solution','Exact solution')
Err = u-uexact;
figure, plot(node,Err,'linewidth',1); legend('Absolute error');

