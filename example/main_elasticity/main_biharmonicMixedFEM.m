clc;clear;close all;
%  -------------- Mesh --------------
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
Nx = 8; Ny = 8; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);

% ------------------ PDE data -------------------
pde = biharmonicdata;

% ----------------- biharmonicMixedFEM ---------------------
maxIt = 4; 
N = zeros(maxIt,1);  h = zeros(maxIt,1);
erruL2 = zeros(maxIt,1);  erruH1 = zeros(maxIt,1);
errwL2 = zeros(maxIt,1);  errwH1 = zeros(maxIt,1);
for k = 1:maxIt
   [node,elem] = uniformrefine(node,elem);     
   bdStruct = setboundary(node,elem);
   [u,w] = biharmonicMixedFEM(node,elem,pde,bdStruct);
   N(k) = size(w,1)+size(u,1);
   h(k) = 1./(sqrt(size(node,1))-1);
   erruL2(k) = getL2error(node,elem,pde.uexact,u);
   erruH1(k) = getH1error(node,elem,pde.Du,u);
   errwL2(k) = getL2error(node,elem,pde.wexact,w);
   errwH1(k) = getH1error(node,elem,pde.Dw,w);
end

% ---------- Plot convergence rates -----------
figure;
subplot(1,2,1);
showrateh2(h,erruH1,1,'-*','|| Du - Du_h||',...
           h,erruL2,1,'k-+','|| u - u_h||');
subplot(1,2,2)
showrateh2(h,errwH1,1,'c-*','|| Dw - Dw_h||',...
           h,errwL2,1,'m-+','|| w - w_h||');    

