clc;clear; close all;
tic;
% -------------- Mesh -----------
a = 0; b = 1;
nel = 10;  N = nel+1; % numbers of elements and nodes
node = linspace(a,b,nel+1)';
elem1D = zeros(nel,2); elem1D(:,1) = 1:N-1; elem1D(:,2) = 2:N;

Neumann = []; Dirichlet = [1,N];
bdStruct = struct('Dirichlet', Dirichlet, 'Neumann', Neumann);

% --------------- PDE ------------
a = -1;  b = 0;  c = 0;
para = struct('a',a, 'b',b, 'c',c);
pde = pde1D(para);

% ----------- FEM1D -----------
u = FEM1D(node,elem1D,pde,bdStruct);

% -------- error analysis -----------
uexact = pde.uexact(node);
figure,plot(node,u,'r-',node,uexact,'k--','linewidth',1);
xlabel('x'); ylabel('u');
legend('Numerical solution','Exact solution')
Err = u-uexact;
figure, plot(node,Err,'linewidth',1); legend('Absolute error');
toc


