clc; clear; close all; 

%% Parameters
maxIt = 5;
N = zeros(maxIt,1);    
h = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1); 
ErrH1 = zeros(maxIt,1);

%% Generate an intitial mesh
[node,elem] = squaremesh([0 1 0 1],0.5);
bdStr = 'x==0';

%% Get PDE data
pde = Poissondatavar;
g_R = @(p) 1 + p(:,1) + p(:,2); % 1 + x + y
pde.g_R = g_R;

%% Finite Element 
i = 3; % 1,2,3
Vh = ['P', num2str(i)];
quadOrder = i+2;
for k = 1:maxIt
    % refine mesh
    [node,elem] = uniformrefine(node,elem);    
    % get the mesh information
    Th = FeMesh2d(node,elem,bdStr);
    % solve the equation
    uh = varPoisson(Th,pde,Vh,quadOrder);
    % record and plot
    N(k) = size(elem,1);
    h(k) = 1/(sqrt(size(node,1))-1);
    if N(k) < 2e3  % show mesh and solution for small size
        figure(1); 
        showresult(node,elem,pde.uexact,uh); 
        drawnow;
    end
    % compute error
    ErrL2(k) = varGetL2Error(Th,pde.uexact,uh,Vh,quadOrder);
    ErrH1(k) = varGetH1Error(Th,pde.Du,uh,Vh,quadOrder);
end

%% Plot convergence rates and display error table
figure(2);
showrateh(h,ErrH1,ErrL2);
fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e');