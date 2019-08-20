clc;clear;close all;

% --------------------- Parameters --------------------
maxN = 1e4;     theta = 0.4;    maxIt = 30;
a1 = 0; b1 = 1; a2 = 0; b2 = 1;

%  -------- Initial mesh and set up PDE data ---------------
Nx = 4; Ny = 4; h1 = 1/Nx; h2 = 1/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);
pde = pdedata();
% ------------- Adaptive Finite Element Method --------------------
for k = 1:maxIt
    % Step 1: SOLVE
    bdFlag = setboundary(node,a1,b1,a2,b2);
    u = Poisson(node,elem,bdFlag,pde);
    figure(1); showmesh(node,elem);
    
    if (size(node,1)>maxN), break; end 
    
    % Step 2: ESTIMATE
    eta = indicator(node,elem,u,pde);
    
    % Step 3: MARK
    elemMarked = mark(elem,eta,theta);
    
    % Step 4: REFINE
    [node,elem] = bisect(node,elem,elemMarked);    
end
% --------- error analysis ---------------
bdFlag = setboundary(node,a1,b1,a2,b2);
u = Poisson(node,elem,bdFlag,pde);
uexact = pde.uexact;
ue = uexact(node);
figure,
subplot(1,2,1), showsolution(node,elem,u); 
zlabel('u');
subplot(1,2,2), showsolution(node,elem,ue); 
zlabel('ue');
Eabs = u-ue;  % Absolute errors
figure,showsolution(node,elem,Eabs); zlim('auto');





