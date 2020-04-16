clc;clear;close all;
%  -------------- Mesh --------------
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
Nx = 4; Ny = 4; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);

% ------------------ PDE data -------------------
pde = biharmonicdata;

% ----------------- biharmonicMixedFEM ---------------------
maxIt = 5; 
N = zeros(maxIt,1);  h = zeros(maxIt,1);
erruL2 = zeros(maxIt,1);  erruH1 = zeros(maxIt,1);
errwL2 = zeros(maxIt,1);  errwH1 = zeros(maxIt,1);

feSpace = 'P1'; 
if strcmpi(feSpace,'P1'), quadOrder = 3; end
if strcmpi(feSpace,'P2'), quadOrder = 4; end
if strcmpi(feSpace,'P3'), quadOrder = 5; end
for k = 1:maxIt
   [node,elem] = uniformrefine(node,elem);     
   bdStruct = setboundary(node,elem);
   Th = getTh(node,elem);
   [u,w] = biharmonicMixedFEM_variational(Th,pde,feSpace,quadOrder);
   NT = size(elem,1);
   h(k) = 1/sqrt(NT);
   erruL2(k) = getL2error(node,elem,u,pde.uexact,feSpace,quadOrder);
   erruH1(k) = getH1error(node,elem,u,pde.Du,feSpace,quadOrder);
   errwL2(k) = getL2error(node,elem,w,pde.wexact,feSpace,quadOrder);
   errwH1(k) = getH1error(node,elem,w,pde.Dw,feSpace,quadOrder);
end

% ---------- Plot convergence rates -----------
figure;
subplot(1,2,1);
showrateh(h, erruL2, erruH1, '|| u - u_h||', '|| Du - Du_h||');
subplot(1,2,2)
showrateh(h, errwL2, errwH1, '|| w - w_h||', '|| Dw - Dw_h||');    

