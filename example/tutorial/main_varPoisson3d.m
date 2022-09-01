clc; clear; close all; 

%% Parameters
maxIt = 5;
N = zeros(maxIt,1);    
h = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1); 
ErrH1 = zeros(maxIt,1);

%% Generate an intitial mesh
%[node,elem3] = cubemesh([0 1 0 1 0 1],0.5);
bdStr = 'x==1';

%% Get the data of the pde
pde = Poissondata3var;
g_R = @(p) 1 + p(:,1) + p(:,2) + p(:,3); % 1 + x + y + z
pde.g_R = g_R;

%% Finite Element Method
i = 1; % 1,2
Vh = ['P', num2str(i)];
quadOrder = i+3;
for k = 1:maxIt
    % refine mesh 
    % WARNING: order of vertices for each cell must be 
    % consistent with the order in iFEM (see line 43 in cubemesh.m).
    %[node,elem3] = uniformrefine3(node,elem3);  
    [node,elem3] = cubemesh([0 1 0 1 0 1],1/(2*k+2));
    % get the mesh information
    Th = FeMesh3d(node,elem3,bdStr);
    Th.solver = 'amg';
    % solve the equation
    uh = varPoisson3d(Th,pde,Vh,quadOrder);
    % record and plot
    N(k) = length(uh);
    h(k) = 1/(size(node,1)^(1/3)-1);
    if N(k) < 2e4  % show mesh and solution for small size
        figure(1);  
        showresult3(node,elem3,pde.uexact,uh);
        drawnow;
    end
    % compute error    
    ErrL2(k) = varGetL2Error3d(Th,pde.uexact,uh,Vh,quadOrder); 
    ErrH1(k) = varGetH1Error3d(Th,pde.Du,uh,Vh,quadOrder);
end

%% Plot convergence rates and display error table
figure(2);
showrateh(h,ErrH1,ErrL2);
fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e');