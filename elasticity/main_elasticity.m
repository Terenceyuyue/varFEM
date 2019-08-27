
clc;clear;close all

% ----------- Mesh and boudary conditions
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
Nx = 10; Ny = 10; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);
bdFlag = setboundary(elem);

% ------------ PDE data ------------
lambda = 1; mu = 1;
para.lambda = lambda; para.mu = mu;
pde = pdedata(para);

% ----------- elasticity --------
u = elasticity(node,elem,pde,bdFlag);

% --------- error analysis -------
u = reshape(u,[],2);
uexact = pde.uexact;  ue = uexact(node);
id = 1;
figure,
subplot(1,2,1), showsolution(node,elem,u(:,id));
zlabel('u');
subplot(1,2,2), showsolution(node,elem,ue(:,id));
zlabel('ue');
Eabs = u-ue;  % Absolute errors
figure,showsolution(node,elem,Eabs(:,id)); zlim('auto');
Err = norm(Eabs)./norm(ue)
