clc; clear; close all;

%% Parameters
maxIt = 5;
N = zeros(maxIt,1);
h = zeros(maxIt,1);
ErruL2 = zeros(maxIt,1);
ErruH1 = zeros(maxIt,1);
ErrwL2 = zeros(maxIt,1);
ErrwH1 = zeros(maxIt,1);

%% Generate an intitial mesh
[node,elem] = squaremesh([0 1 0 1],0.25);

%% Get the data of the pde
pde = biharmonicdatavar;  % slight difference from biharmonicdata

%% Finite Element Method
i = 1; % 1,2,3
Vh = ['P', num2str(i)];
quadOrder = i+2;
Vhvec = repmat( {Vh}, 1, 2 ); % v = [v1,v2]
for k = 1:maxIt
    % refine mesh
    [node,elem] = uniformrefine(node,elem);
    % get the mesh information
    Th = FeMesh2d(node,elem);
    % solve the equation
    [uh,wh] = varBiharmonicMixedFEM(Th,pde,Vhvec,quadOrder);
    % record and plot
    N(k) = size(elem,1);
    h(k) = 1/(sqrt(size(node,1))-1);
    if N(k) < 2e3  % show mesh and solution for small size
        figure(1);  
        showresult(node,elem,pde.uexact,uh);
        drawnow;
    end
    % compute error
    ErruL2(k) = varGetL2Error(Th,pde.uexact,uh,Vh,quadOrder);
    ErruH1(k) = varGetH1Error(Th,pde.Du,uh,Vh,quadOrder);
    ErrwL2(k) = varGetL2Error(Th,pde.wexact,wh,Vh,quadOrder);
    ErrwH1(k) = varGetH1Error(Th,pde.Dw,wh,Vh,quadOrder);
end

%% Plot convergence rates and display error table
figure(2);
subplot(1,2,1)
showrateErr(h,ErruH1,ErruL2);
subplot(1,2,2)
showrateErr(h,ErrwH1,ErrwL2);
fprintf('\n');
disp('Table: Error of uh')
colname = {'#Dof','h','||u-u_h||','|Du-Du_h|_1'};
disptable(colname,N,[],h,'%0.3e',ErruL2,'%0.5e',ErruH1,'%0.5e');

fprintf('\n');
disp('Table: Error of wh')
colname = {'#Dof','h','||w-w_h||','||Dw-Dw_h||'};
disptable(colname,N,[],h,'%0.3e',ErrwL2,'%0.5e',ErrwH1,'%0.5e');