clc; clear; close all;

%% Parameters
maxIt = 5;
N = zeros(maxIt,1);
h = zeros(maxIt,1);
ErruL2 = zeros(maxIt,1);
ErruH1 = zeros(maxIt,1);
ErrpL2 = zeros(maxIt,1);
ErrpH1 = zeros(maxIt,1);

%% Generate an intitial mesh
[node,elem] = squaremesh([0 1 0 1],0.5);

%% Get the data of the pde
nu = 10;
pde = NavierStokesdata(nu);

%% Finite Element Method
Vh = {'P2','P2','P1'}; % v = [v1,v2,q] --> [v1,v2,v3]
quadOrder = 5;
for k = 1:maxIt
    % refine mesh
    [node,elem] = uniformrefine(node,elem);
    % get the mesh information
    Th = FeMesh2d(node,elem); 
    % solve the equation
    [uh,ph] = varNavierStokes(Th,pde,Vh,quadOrder);
    % record and plot
    N(k) = size(elem,1);
    h(k) = 1/(sqrt(size(node,1))-1);
    if N(k) < 2e3  % show mesh and solution for small size
        figure(1);  
        showresult(node,elem,pde.uexact,uh);
        drawnow;
    end
    % compute error
    tru = eye(2); trDu = eye(4);
    erruL2 = zeros(1,2);  erruH1 = zeros(1,2); % square
    for id = 1:2
        uid = uh(:,id);
        u = @(pz) pde.uexact(pz)*tru(:, id);
        Du = @(pz) pde.Du(pz)*trDu(:, 2*id-1:2*id);
        erruL2(:,id) = varGetL2Error(Th,u,uid,Vh{1},quadOrder);
        erruH1(:,id) = varGetH1Error(Th,Du,uid,Vh{1},quadOrder);
    end
    
    ErruL2(k) = sqrt(sum(erruL2.^2,2));
    ErruH1(k) = sqrt(sum(erruH1.^2,2));
    
    ErrpL2(k) = varGetL2Error(Th,pde.pexact,ph,Vh{3},quadOrder);
    ErrpH1(k) = varGetH1Error(Th,pde.Dpexact,ph,Vh{3},quadOrder);
end

%% Plot convergence rates and display error table
figure;
pe = pde.pexact(node);
x = 0:0.01:1; y = x;
subplot(1,2,1),
varcontourf(x,y,node,elem,ph,20);
title('Numerical pressure');
subplot(1,2,2),
varcontourf(x,y,node,elem,pe,20);
title('Exact pressure');

figure;
subplot(1,2,1)
showrateErr(h,ErruH1,ErruL2);
subplot(1,2,2)
showrateErr(h,ErrpH1,ErrpL2);
fprintf('\n');
disp('Table: Error of uh')
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||', '||p-p_h||'};
disptable(colname,N,[],h,'%0.3e',ErruL2,'%0.5e',ErruH1,'%0.5e',ErrpL2,'%0.5e');

fprintf('\n');
disp('Table: Error of ph')
colname = {'#Dof','h','||p-p_h||','||Dp-Dp_h||'};
disptable(colname,N,[],h,'%0.3e',ErrpL2,'%0.5e',ErrpH1,'%0.5e');