clc; clear; close all;
%varFEM1D solves the two-point boundary value problem
%
%     -au'' + bu' + cu = f  in \Omega = (a,b), with
%     Dirichlet boundary condition u = g_D  on \Gamma_D = {a}, {b} or {a,b},
%     Neumann boundary condition   u = u' on \Gamma_N = {a,b} - \Gamma_D
%

%% Parameters
maxIt = 5;
N = zeros(maxIt,1);
h = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);

%% Generate an intitial mesh
a = 0; b = 1;
nel = 4;  % numbers of elements
node = linspace(a,b,nel+1)';
elem1d = [(1:nel)', (2:nel+1)'];
bdNeumann = 'x==0';

%% Get the data of the pde
pde = pdedata1Dvar;

%% Finite Element Method
i = 1; % 1,2,3
Vh = ['P', num2str(i)];
quadOrder = i+3;
for k = 1:maxIt
    % step 1: refine mesh 
    [node,elem1d] = uniformrefine1(node,elem1d);
    
    % step 2: get the mesh information
    Th = FeMesh1d(node,elem1d,bdNeumann);
    
    % step 3: assemble stiffness matrix
    Coef  = {pde.a,   pde.b,   pde.c};
    Test  = {'v.dx', 'v.val',  'v.val'};
    Trial = {'u.dx',  'u.dx',  'u.val'};
    kk = assem1d(Th,Coef,Test,Trial,Vh,quadOrder); 
    
    % step 4: assemble the right hand side
    Coef = pde.f; Test = 'v.val';
    ff = assem1d(Th,Coef,Test,[],Vh,quadOrder);
    
    % step 5: apply boundary conditions
    uh = apply1d(Th,kk,ff,pde);
    
    % step 6: show solution
    N(k) = size(elem1d,1); h(k) = 1/size(elem1d,1);
    if N(k) < 1e2  % show solution for small size
        figure(1);
        [node1,id] = sort(node);
        plot(node1,uh(id),'r-',node1,pde.uexact(node1), ...
            'k*','linewidth',2);
        pause(0.5);
    end
    
    % step 7: compute errors
    ErrL2(k) = varGetL2Error1d(Th,pde.uexact,uh,Vh,quadOrder);
    ErrH1(k) = varGetH1Error1d(Th,pde.Du,uh,Vh,quadOrder);
end

%% Plot convergence rates and display error table
figure(2);
showrateh(h,ErrH1,ErrL2);
fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e');