clc; close all;
clear variables;

%% Parameters
maxIt = 5;
N = zeros(maxIt,1);
h = zeros(maxIt,1);
erruL2 = zeros(maxIt,1);
erruH1 = zeros(maxIt,1);
errwL2 = zeros(maxIt,1);
errwH1 = zeros(maxIt,1);

%% Generate an intitial mesh
[node,elem] = squaremesh([0 1 0 1],0.25);

%% Get the data of the pde
pde = biharmonicdatavar;

%% Finite Element Method
for k = 1:maxIt
    % refine mesh
    [node,elem] = uniformrefine(node,elem);
    % get the mesh information
    Th = getTh2D(node,elem);
    % solve the equation
    [uh,wh] = varBiharmonicMixedFEM_block(Th,pde);
    % record and plot
    N(k) = size(elem,1);
    h(k) = 1/(sqrt(size(node,1))-1);
    if N(k) < 2e3  % show mesh and solution for small size
        figure(1);  
        showresult(node,elem,pde.exactu,uh);
        drawnow;
    end
    erruL2(k) = getL2error(node,elem,pde.exactu,uh);
    erruH1(k) = getH1error(node,elem,pde.Du,uh);
    errwL2(k) = getL2error(node,elem,pde.exactw,wh);
    errwH1(k) = getH1error(node,elem,pde.Dw,wh);
end

%% Plot convergence rates and display error table
figure(2);
subplot(1,2,1)
showrateh(h,erruH1,erruL2);
subplot(1,2,2)
showrateh(h,errwH1,errwL2);

fprintf('\n');
disp('Table: Error of uh')
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||'};
disptable(colname,N,[],h,'%0.3e',erruL2,'%0.5e',erruH1,'%0.5e');

fprintf('\n');
disp('Table: Error of wh')
colname = {'#Dof','h','||w-w_h||','||Dw-Dw_h||'};
disptable(colname,N,[],h,'%0.3e',errwL2,'%0.5e',errwH1,'%0.5e');

