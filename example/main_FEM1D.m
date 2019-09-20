clc;clear; close all;
tic;
% -------------- Mesh -----------
a = 0; b = 1;
nel = 10;  N = nel+1; % numbers of elements and nodes
node = linspace(a,b,nel+1)';
elem = zeros(nel,2); elem(:,1) = 1:N-1; elem(:,2) = 2:N;

Neumann = []; Dirichlet = [1,N];
bdStruct = struct('Dirichlet', Dirichlet, 'Neumann', Neumann);

% --------------- PDE ------------
acoef = -1;  bcoef = 1;  ccoef = 2;
para = struct('acoef',acoef, 'bcoef',bcoef, 'ccoef',ccoef);
pde = pdedata1D(para);

% ----------- elasticity --------
u = FEM1D(node,elem,pde,bdStruct);

% -------- error analysis -----------
uexact = pde.uexact(node);
figure,plot(node,u,'r-',node,uexact,'k--','linewidth',1);
xlabel('x'); ylabel('u');
legend('Numerical solution','Exact solution')
Err = u-uexact;
figure, plot(node,Err,'linewidth',1); legend('Absolute error');
toc


