clc; close all;
clear variables;

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
elem = [(1:nel)', (2:nel+1)'];
bdNeumann = 'abs(x-1)<1e-4';

%% Get the data of the pde
pde = pdedata1Dvar;

%% Finite Element Method
Vh = 'P3';
if strcmpi(Vh,'P1'), quadOrder = 4; end
if strcmpi(Vh,'P2'), quadOrder = 5; end
if strcmpi(Vh,'P3'), quadOrder = 6; end
for k = 1:maxIt
    % refine mesh
    [node,elem] = uniformrefine1(node,elem);
    % get the mesh information
    Th = getTh1D(node,elem,bdNeumann);
    % solve the equation
    uh = varFEM1D(Th,pde,Vh,quadOrder);
    % record
    N(k) = size(elem,1);
    h(k) = 1/size(elem,1);
    if N(k) < 1e2  % show solution for small size
        figure(1);
        [node1,id] = sort(node);
        plot(node1,uh(id),'r-',node1,pde.uexact(node1), ...
            'k*','linewidth',2);
        pause(1);
    end
    % compute error
    ErrL2(k) = getL2error1(node,elem,pde.uexact,uh,quadOrder);
    ErrH1(k) = getH1error1(node,elem,pde.Du,uh,quadOrder);
end

%% Plot convergence rates and display error table
figure(2);
showrateh(h,ErrH1,ErrL2);
fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e');

