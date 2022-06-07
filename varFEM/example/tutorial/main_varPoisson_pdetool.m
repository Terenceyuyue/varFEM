clc; clear; close all;

%% Parameters
maxIt = 4;
N = zeros(maxIt,1);
h = zeros(maxIt,1);
errL2 = zeros(maxIt,1);
errH1 = zeros(maxIt,1);

%% Generate an intitial mesh
g = 'circleg'; %'lshapeg';
[p,e,t] = initmesh(g,'hmax',0.5);
node = p'; elem = t(1:3,:)';
bdNeumann = 'x>0'; % string for Neumann

%% Get PDE data
pde = Poissondatavar;
g_R = @(p) 1 + p(:,1) + p(:,2); % 1 + x + y
pde.g_R = g_R;

%% Finite Element Method
i = 1; % 1,2,3
Vh = ['P', num2str(i)];
quadOrder = i+2;
for k = 1:maxIt
    % refine mesh
    [p,e,t] = refinemesh(g,p,e,t);
    % get the mesh information
    node = p'; elem = t(1:3,:)';
    Th = FeMesh2d(node,elem,bdNeumann);
    % solve the equation
    uh = varPoisson(Th,pde,Vh,quadOrder);
    % record and plot
    N(k) = size(elem,1);
    h(k) = 1/(sqrt(size(node,1))-1);
    if N(k) < 2e3  % show mesh and solution for small size
        figure(1);  
        showresult(node,elem,pde.uexact,uh);
        pause(1);
    end
    % compute error
    errL2(k) = getL2error(node,elem,pde.uexact,uh);
    errH1(k) = getH1error(node,elem,pde.Du,uh);
end

%% Plot convergence rates and display error table
figure(2);
showrateh(h,errH1,errL2);
fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||'};
disptable(colname,N,[],h,'%0.3e',errL2,'%0.5e',errH1,'%0.5e');