clc;clear; close all;
% ----------- Mesh and boundary conditions -----------
a = 0; b = 1;
nel = 4;  N = nel+1; % numbers of elements and nodes
node = linspace(a,b,nel+1)';
elem1D = zeros(nel,2); elem1D(:,1) = 1:N-1; elem1D(:,2) = 2:N;

bdNeumann = 'abs(x-1)<=1e-4';

% --------------------- PDE ---------------------------
pde = pdedata1D;

% ------------------------ Poisson ------------------------
maxIt = 5;
N = zeros(maxIt,1);    h = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1); ErrH1 = zeros(maxIt,1);

feSpace = 'P3';
if strcmpi(feSpace,'P1'), quadOrder = 4; end
if strcmpi(feSpace,'P2'), quadOrder = 5; end
if strcmpi(feSpace,'P3'), quadOrder = 6; end
for k = 1:maxIt
    [node,elem1D] = uniformrefine1D(node,elem1D);  
    bdStruct = setboundary1D(node,elem1D,bdNeumann);
    Th.node = node; Th.elem1D = elem1D; Th.bdStruct = bdStruct;
    uh = FEM1D_variational(Th,pde,feSpace,quadOrder); 
    nel = size(elem1D,1);
    h(k) = 1/nel;
    ErrL2(k) = getL2error1D(node,elem1D,uh,pde.uexact,feSpace,quadOrder);
    ErrH1(k) = getH1error1D(node,elem1D,uh,pde.Du,feSpace,quadOrder);
end

% -------------------- Show rate -----------------------
figure,showrateh(h,ErrL2,ErrH1);

% [node,id] = sort(node);
% uh = uh(id);
% figure,plot(node,uh,'k-',node,pde.uexact(node),'r*','linewidth',1);
% xlabel('x'); ylabel('u');
% legend('Numerical solution','Exact solution')
