clc;clear;close all;

%% Parameters
maxN = 1e4;     theta = 0.4;    maxIt = 100;
h = zeros(maxIt,1);  N = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);
etaN = zeros(maxIt,1);

%% Generate an initial mesh
[node,elem] = squaremesh([0 1 0 1],1/5);
bdStr = [];

%% Get PDE data
pde = Poissondata_afem();

%% Adaptive Finite Element Method 
Vh = 'P2'; quadOrder = 7;
for k = 1:maxIt
    % Step 1: SOLVE
    % get the mesh information
    Th = FeMesh2d(node,elem,bdStr);
    % solve the equation
    uh = varPoisson(Th,pde,Vh,quadOrder);
    figure(1); 
    showresult(node,elem,pde.uexact,uh);
    drawnow;
    ErrL2(k) = varGetL2Error(Th,pde.uexact,uh,Vh,quadOrder);
    ErrH1(k) = varGetH1Error(Th,pde.Du,uh,Vh,quadOrder);
    
    % Step 2: ESTIMATE
    eta = Poisson_indicator(Th,uh,pde,Vh,quadOrder);
    
    % Step 3: MARK
    elemMarked = mark(elem,eta,theta);
    
    % Step 4: REFINE
    [node,elem] = bisect(node,elem,elemMarked);

    % Record
    etaN(k) = norm(eta);
    N(k) = length(uh);  
    h(k) = 1/sqrt(size(elem,1));
    
    if (size(node,1)>maxN) || (k==maxIt)
        uh = varPoisson(Th,pde,Vh,quadOrder);
        step = k
        break;
    end    
end

figure,
plot((1:step),etaN(1:step),'k.-','linewidth',1);
xlabel('k'); ylabel('\eta (u_h)');

figure,
id = 15;
h = 1./sqrt(N(id:step));
showratehh(h,etaN(id:step),'r-*','\eta (u_h)', ErrH1(id:step), 'b-s','|u-u_h|_1')


