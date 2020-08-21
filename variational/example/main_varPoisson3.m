clc; close all; 
clear variables;

%% Parameters
maxIt = 4;
N = zeros(maxIt,1);    
h = zeros(maxIt,1);
errL2 = zeros(maxIt,1); 
errH1 = zeros(maxIt,1);

%% Generate an intitial mesh
[node,elem] = cubemesh([0 1 0 1 0 1],1);
bdNeumann = 'abs(x-1)<1e-4';

%% Get the data of the pde
pde = Poissondata3var;
g_R = @(p) 1 + p(:,1) + p(:,2) + p(:,3); % 1 + x + y + z
pde.g_R = g_R;

%% Finite Element Method
Vh = 'P2';
if strcmpi(Vh,'P1'), quadOrder = 3; end
if strcmpi(Vh,'P2'), quadOrder = 4; end
for k = 1:maxIt
    % refine mesh 
    % WARNING: order of vertices for each cell must be 
    % consistent with the order in iFEM (see line 43 in cubemesh.m).
    [node,elem] = uniformrefine3(node,elem);  
    % get the mesh information
    Th = getTh3D(node,elem,bdNeumann);
    % solve the equation
    uh = varPoisson3(Th,pde,Vh,quadOrder);
    % record and plot
    N(k) = length(uh);
    h(k) = 1/(size(node,1)^(1/3)-1);
    if N(k) < 2e4  % show mesh and solution for small size
        figure(1);  
        showresult3(node,elem,pde.exactu,uh);
        pause(0.1);
    end
    % compute error    
    errL2(k) = getL2error3(node,Th.elem,pde.exactu,uh); 
    errH1(k) = getH1error3(node,Th.elem,pde.Du,uh);
end

%% Plot convergence rates and display error table
figure(2);
showrateh(h,errH1,errL2);
fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1'};
disptable(colname,N,[],h,'%0.3e',errL2,'%0.5e',errH1,'%0.5e');


