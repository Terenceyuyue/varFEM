clc; clear; close all; 

%% Parameters
maxIt = 3;
N = zeros(maxIt,1);    
h = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1); 
ErrH1 = zeros(maxIt,1);

%% Generate an intitial mesh
[node,elem3] = cubemesh([0 1 0 1 0 1],0.25);
bdStr = 'x==1'; % string for Neumann

%% Get the data of the pde
lambda = 1; mu = 1;
para.lambda = lambda; para.mu = mu;
pde = elasticitydata3d(para);

%% Finite Element Method
i = 1; % 1,2
Vh = ['P', num2str(i)];
quadOrder = i+3;
Vhvec = repmat( {Vh}, 1, 3 ); % v = [v1,v2,v3]
for k = 1:maxIt
    % refine mesh 
    % WARNING: order of vertices for each cell must be 
    % consistent with the order in iFEM (see line 43 in cubemesh.m).
    [node,elem3] = uniformrefine3(node,elem3); 
    % get the mesh information
    Th = FeMesh3d(node,elem3,bdStr); 
    Th.solver = 'amg';
    % solve the equation
    uh = varElasticity3d(Th,pde,Vhvec,quadOrder);
    uh = reshape(uh,[],3);
    % record and plot
    N(k) = 3*size(uh,1);
    h(k) = 1/(size(node,1)^(1/3)-1);
    if N(k) < 2e4  % show mesh and solution for small size
        figure(1);  
        tru = eye(3); i = 2;
        ue = @(pz) pde.uexact(pz)*tru(:,i);
        showresult3(node,elem3,ue,uh(:,i));
        drawnow;
    end
    % compute error
    tru = eye(3); trDu = eye(3*3);
    errL2 = zeros(1,3);  errH1 = zeros(1,3); % square
    for id = 1:3
        uhid = uh(:,id);
        u = @(pz) pde.uexact(pz)*tru(:, id);
        Du = @(pz) pde.Du(pz)*trDu(:, 3*id-2:3*id);
        errL2(id) = varGetL2Error3d(Th,u,uhid,Vh,quadOrder);
        errH1(id) = varGetH1Error3d(Th,Du,uhid,Vh,quadOrder);
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