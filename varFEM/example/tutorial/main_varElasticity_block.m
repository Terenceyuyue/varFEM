clc; clear; close all; 

%% Parameters
maxIt = 4;
N = zeros(maxIt,1);    
h = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1); 
ErrH1 = zeros(maxIt,1);

%% Generate an intitial mesh
[node,elem] = squaremesh([0 1 0 1],0.25);
bdNeumann = 'y==0 | x==1'; % string for Neumann

%% Get the data of the pde
lambda = 1; mu = 1;
para.lambda = lambda; para.mu = mu;
pde = elasticitydatavar(para);

%% Finite Element Method
for k = 1:maxIt
    % refine mesh
    [node,elem] = uniformrefine(node,elem);
    % get the mesh information
    Th = FeMesh2d(node,elem,bdNeumann); 
    % solve the equation
    uh = varElasticity_block(Th,pde);
    uh = reshape(uh,[],2);
    % record and plot
    N(k) = size(elem,1);
    h(k) = 1/(sqrt(size(node,1))-1);
    if N(k) < 2e3  % show mesh and solution for small size
        figure(1);  
        showresult(node,elem,pde.uexact,uh(:,1));
        pause(1);
    end
    % compute error
    tru = eye(2); trDu = eye(4);
    errL2 = zeros(1,2);  errH1 = zeros(1,2); % square
    for id = 1:2
        uid = uh(:,id);
        u = @(pz) pde.uexact(pz)*tru(:, id);
        Du = @(pz) pde.Du(pz)*trDu(:, 2*id-1:2*id);
        Vh = 'P1';  quadOrder = 5;
        errL2(:,id) = varGetL2Error(Th,u,uid,Vh,quadOrder);
        errH1(:,id) = varGetH1Error(Th,Du,uid,Vh,quadOrder);
    end
    
    ErrL2(k) = sqrt(sum(errL2.^2,2));
    ErrH1(k) = sqrt(sum(errH1.^2,2));
end

%% Plot convergence rates and display error table
figure(2);
showrateh(h,ErrH1,ErrL2);
fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e');
