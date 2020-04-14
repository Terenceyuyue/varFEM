clc;clear;close all;
% --------------- Mesh and boundary conditions ------------------
g = 'circleg'; %'lshapeg';
[p,e,t] = initmesh(g,'hmax',0.5);
node = p'; elem = t(1:3,:)';

bdNeumann = 'x>0'; % string for Neumann

% ------------------------- PDE data ------------------------
pde = Poissondata2();
g_R = @(p) 1 + p(:,1) + p(:,2); % 1 + x + y
pde.g_R = g_R;

% ------------------------ Poisson ------------------------
maxIt = 5;
N = zeros(maxIt,1);    h = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1); ErrH1 = zeros(maxIt,1);

feSpace = 'P2';
if strcmpi(feSpace,'P1'), quadOrder = 3; end
if strcmpi(feSpace,'P2'), quadOrder = 4; end
if strcmpi(feSpace,'P3'), quadOrder = 5; end
for k = 1:maxIt
    [p,e,t] = refinemesh(g,p,e,t);
    node = p'; elem = t(1:3,:)';
    %figure(1), showmesh(node,elem), hold on, pause(0.6)
    bdStruct = setboundary(node,elem,bdNeumann);
    Th.node = node; Th.elem = elem; Th.bdStruct = bdStruct;    
    uh = Poisson_variational(Th,pde,feSpace,quadOrder);
    NT = size(elem,1);
    h(k) = 1/sqrt(NT);
    ErrL2(k) = getL2error(node,elem,uh,pde.uexact,feSpace,quadOrder);
    ErrH1(k) = getH1error(node,elem,uh,pde.Du,feSpace,quadOrder);
end

% -------------------- Show rate -----------------------
figure,showrateh(h,ErrL2,ErrH1);