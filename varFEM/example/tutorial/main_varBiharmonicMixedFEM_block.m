clc; clear; close all;

%% Parameters
maxIt = 4;
N = zeros(maxIt,1);
h = zeros(maxIt,1);
ErruL2 = zeros(maxIt,1);
ErruH1 = zeros(maxIt,1);
ErrwL2 = zeros(maxIt,1);
ErrwH1 = zeros(maxIt,1);

%% Generate an intitial mesh
[node,elem] = squaremesh([0 1 0 1],0.25);

%% Get the data of the pde
pde = biharmonicdatavar;

%% Finite Element Method
for k = 1:maxIt
    % refine mesh
    [node,elem] = uniformrefine(node,elem);
    % get the mesh information
    Th = FeMesh2d(node,elem);
    % solve the equation
    [uh,wh] = varBiharmonicMixedFEM_block(Th,pde);
    % record and plot
    N(k) = size(elem,1);
    h(k) = 1/(sqrt(size(node,1))-1);
    if N(k) < 2e3  % show mesh and solution for small size
        figure(1);  
        showresult(node,elem,pde.uexact,uh);
        drawnow;
    end
    ErruL2(k) = getL2error(node,elem,pde.uexact,uh);
    ErruH1(k) = getH1error(node,elem,pde.Du,uh);
    ErrwL2(k) = getL2error(node,elem,pde.wexact,wh);
    ErrwH1(k) = getH1error(node,elem,pde.Dw,wh);
end

%% Plot convergence rates and display error table
figure(2);
subplot(1,2,1)
showrateh(h,ErruH1,ErruL2);
subplot(1,2,2)
showrateh(h,ErrwH1,ErrwL2);

fprintf('\n');
disp('Table: Error of uh')
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||'};
disptable(colname,N,[],h,'%0.3e',ErruL2,'%0.5e',ErruH1,'%0.5e');

fprintf('\n');
disp('Table: Error of wh')
colname = {'#Dof','h','||w-w_h||','||Dw-Dw_h||'};
disptable(colname,N,[],h,'%0.3e',ErrwL2,'%0.5e',ErrwH1,'%0.5e');