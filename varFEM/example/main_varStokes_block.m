clc; close all;
clear variables;

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
pde = Stokesdatavar;

%% Finite Element Method

Vh = {'P2','P2','P1'}; quadOrder = 4;
for k = 1:maxIt
    % refine mesh
    [node,elem] = uniformrefine(node,elem);
    % get the mesh information
    Th = getTh2D(node,elem);
    % solve the equation
    [uh,ph] = varStokes_block(Th,pde,Vh,quadOrder);
    % record and plot
    N(k) = size(node,1);
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
        erruL2(:,id) = getL2error(node,elem,u,uid,quadOrder);
        erruH1(:,id) = getH1error(node,elem,Du,uid,quadOrder);
    end
    
    ErruL2(k) = sqrt(sum(erruL2.^2,2));
    ErruH1(k) = sqrt(sum(erruH1.^2,2));
    
    ErrpL2(k) = getL2error(node,elem,pde.pexact,ph,quadOrder);
    ErrpH1(k) = getH1error(node,elem,pde.Dpexact,ph,quadOrder);
end

%% Plot convergence rates and display error table
figure(2);
subplot(1,2,1)
showrateh(h,ErruH1,ErruL2);
subplot(1,2,2)
showrateh(h,ErrpH1,ErrpL2);
fprintf('\n');
disp('Table: Error of uh')
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||'};
disptable(colname,N,[],h,'%0.3e',ErruL2,'%0.5e',ErruH1,'%0.5e');

fprintf('\n');
disp('Table: Error of ph')
colname = {'#Dof','h','||p-p_h||','||Dp-Dp_h||'};
disptable(colname,N,[],h,'%0.3e',ErrpL2,'%0.5e',ErrpH1,'%0.5e');
