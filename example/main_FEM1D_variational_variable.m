clc;clear; close all;
tic;
% -------------- Mesh and boundary conditions -----------
a = 0; b = 1;
nel = 10;  N = nel+1; % numbers of elements and nodes
node = linspace(a,b,nel+1)';
elem = zeros(nel,2); elem(:,1) = 1:N-1; elem(:,2) = 2:N;
Th.node = node; Th.elem = elem;

Neumann = []; Dirichlet = [1,N];
bdStruct = struct('Dirichlet', Dirichlet, 'Neumann', Neumann);

% --------------- PDE ------------
acoef = @(x) 10+x;  bcoef = @(x) 1+x;  ccoef = 2;
para = struct('acoef',acoef, 'bcoef',bcoef, 'ccoef',ccoef);
pde = pdedata1D_variable(para);

% --------------------- Assemble the matrix ----------------
acoef1 = @(x) -acoef(x);
Coef  = {acoef1,   bcoef,     ccoef};
Trial = {'u.dx',  'u.dx',    'u.val'};
Test  = {'v.dx', 'v.val',    'v.val'};
kk = Bilinear1D(Th,Coef,Trial,Test);

% --------- Assemble the vector ------------
ff = Linear1D(Th,pde.f);

% --------- Apply boundary conditions ---------
u = Applyboundary1D(Th,kk,ff,acoef,pde,bdStruct);

% -------- error analysis -----------
uexact = pde.uexact(node);
figure,plot(node,u,'r-',node,uexact,'k--','linewidth',1);
xlabel('x'); ylabel('u');
legend('Numerical solution','Exact solution')
Err = u-uexact;
figure, plot(node,Err,'linewidth',1); legend('Absolute error');
toc
