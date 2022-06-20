clc; clear; close all; 

%% Parameters
maxIt = 3;
N = zeros(maxIt,1);    
h = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1); 
ErrH1 = zeros(maxIt,1);

%% Generate an intitial mesh
[node,elem] = squaremesh([0 1 0 1],0.5);
bdNeumann = 'x==0';

%% Get PDE data
pde = heatData();
pde.t0 = 0; pde.tf = 1; 
pde.tn = 0.05; % compute values at t = tn

%% Finite Element 
iq = 1; % 1,2,3
Vh = ['P', num2str(iq)];
quadOrder = iq+2;
for k = 1:maxIt
    % refine mesh
    [node,elem] = uniformrefine(node,elem);    
    % get the mesh information
    Th = FeMesh2d(node,elem,bdNeumann);
    pde.Nt = (2^(1+k))^(iq+1);  % tau = c*h^(iq+1)
    % solve the equation
    [uh,ErrH1k,ErrL2k] = varHeat(Th,pde,Vh,quadOrder);
    % record and plot
    N(k) = size(elem,1);
    h(k) = 1/(sqrt(size(node,1))-1);
    if N(k) < 2e3  % show mesh and solution for small size
        figure(1); 
        ue = @(p) pde.uexact(p,pde.tn);
        showresult(node,elem,ue,uh); 
        drawnow;
    end
    % compute error
    ErrL2(k) = ErrL2k;
    ErrH1(k) = ErrH1k;
end

%% Plot convergence rates and display error table
figure(2);
showrateh(h,ErrH1,ErrL2);
fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e');