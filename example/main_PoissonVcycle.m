clc;clear;close all;
%  -------------- Initial mesh --------------
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
Nx = 5; Ny = 5; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);
J = 3; % solve it directly when J<=1

% ------ Final mesh, prolongation matrices and boundary conditions -------
[node,elem,Pro,Res] = MeshPoissonVcycle(node,elem,J);

bdNeumann = 'abs(x-1)<1e-4'; % string for Neumann
bdStruct = setboundary(node,elem,bdNeumann);

% ------------------ PDE data -------------------
% pde = struct('uexact',@uexact, 'f',@f, 'g_N',@g_N,  'g_D',@g_D);
pde = Poissondata();

% ------------------ Poisson ---------------------
u = PoissonVcycle(node,elem,pde,bdStruct,Pro,Res);

% --------- error analysis ---------------
uexact = pde.uexact;
ue = uexact(node);
figure,
subplot(1,2,1), showsolution(node,elem,u); 
subplot(1,2,2), showsolution(node,elem,ue);
Eabs = u-ue;  % Absolute errors
figure,showsolution(node,elem,Eabs); zlim('auto');
format shorte
Err = norm(u-ue)/norm(ue)