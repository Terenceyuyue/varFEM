clc;clear;close all;
%  -------------- Mesh and boundary conditions --------------
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
Nx = 2; Ny = 2; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);

bdNeumann = 'abs(x-1)<1e-4'; % string for Neumann

% ----------------------- PDE data ------------------------
pde = Poissondata();

% ------------------------ Poisson ------------------------
maxIt = 5;
N = zeros(maxIt,1);    h = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1); ErrH1 = zeros(maxIt,1);

for k = 1:maxIt
    [node,elem] = uniformrefine(node,elem);
    bdStruct = setboundary(node,elem,bdNeumann);
    uh = Poisson(node,elem,pde,bdStruct);
    NT = size(elem,1);
    h(k) = 1/sqrt(NT);
    ErrL2(k) = getL2error(node,elem,uh,pde.uexact);
    ErrH1(k) = getH1error(node,elem,uh,pde.Du);
end

% -------------------- Show rate -----------------------
showrateh(h,ErrL2,ErrH1);


