
clc;clear;close all;
%  -------------- Mesh and PDE --------------
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
Nx = 50; Ny = 50; h1 = 1/Nx; h2 = 1/Ny;
[node,elem] = squaremesh([0 1 0 1],h1,h2);
pde = pdedata();
bdFlag = setboundary(node,a1,b1,a2,b2);

% --------------- solve Poisson problem -------
u = Poisson(node,elem,bdFlag,pde);

% --------- error analysis ---------------
uexact = pde.uexact;
ue = uexact(node);
figure,
subplot(1,2,1), showsolution(node,elem,u); 
zlabel('u');
subplot(1,2,2), showsolution(node,elem,ue); 
zlabel('ue');
Eabs = u-ue;  % Absolute errors
figure,showsolution(node,elem,Eabs); zlim('auto');


