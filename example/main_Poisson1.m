clc;clear;close all;
%  -------------- Mesh and boundary conditions --------------
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
Nx = 5; Ny = 5; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);

bdNeumann = 'abs(x-1)<1e-4'; % string for Neumann

% ------------------ PDE data -------------------
% pde = struct('uexact',@uexact, 'f',@f, 'g_N',@g_N,  'g_D',@g_D);
pde = Poissondata();

% ------------------ Poisson ---------------------
maxIt = 5; 
N = zeros(maxIt,1);    h = zeros(maxIt,1);
erruL2 = zeros(maxIt,1); erruH1 = zeros(maxIt,1);

for k = 1:maxIt
   [node,elem] = uniformrefine(node,elem);  
   bdStruct = setboundary(node,elem,bdNeumann);
   u = Poisson(node,elem,pde,bdStruct);
   N(k) = size(u,1);
   h(k) = 1./(sqrt(size(node,1))-1);
   erruL2(k) = getL2error(node,elem,pde.uexact,u);
   erruH1(k) = getH1error(node,elem,pde.Du,u);
end

showrateh2(h,erruH1,1,'-*','|| Du - Du_h||',...
           h,erruL2,1,'k-+','|| u - u_h||');
