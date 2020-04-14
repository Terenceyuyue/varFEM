clc;clear;close all;
%  -------------- Mesh and boundary conditions --------------
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
Nx = 10; Ny = 10; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);

bdNeumann = 'abs(x-1)<1e-4'; % string for Neumann
bdStruct = setboundary(node,elem,bdNeumann);

% ------------------ PDE data -------------------
pde = Poissondata();

% ------------------ Poisson ---------------------
uh = Poisson(node,elem,pde,bdStruct);

% --------- error analysis ---------------
uexact = pde.uexact;
ue = uexact(node);
figure,
subplot(1,2,1), showsolution(node,elem,uh); 
subplot(1,2,2), showsolution(node,elem,ue);

% L2 and H1 errors
format shorte
ErrL2 = getL2error(node,elem,uh,uexact)
ErrH1 = getH1error(node,elem,uh,pde.Du)