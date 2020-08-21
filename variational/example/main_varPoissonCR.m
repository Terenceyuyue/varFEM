clc; close all;
clear variables;

%% Parameters
maxIt = 5;
N = zeros(maxIt,1);
h = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);

%% Generate an intitial mesh
[node,elem] = squaremesh([0 1 0 1],0.5);
bdNeumann = 'abs(x-0)<1e-4';

%% Get the data of the pde
pde = Poissondatavar;
g_R = @(p) 1 + p(:,1) + p(:,2); % 1 + x + y
pde.g_R = g_R;

%% Finite Element Method
Vh = 'CR';
if strcmpi(Vh,'CR'), quadOrder = 3; end
for k = 1:maxIt
    % refine mesh
    [node,elem] = uniformrefine(node,elem);
    % get the mesh information
    Th = getTh2D(node,elem,bdNeumann);
    % solve the equation
    uh = varPoissonCR(Th,pde,Vh,quadOrder);
    % record and plot
    N(k) = size(elem,1);
    h(k) = 1/(sqrt(size(node,1))-1);
    if N(k) < 2e3  % show mesh and solution for small size
        figure(1);
        NT = size(elem,1);        
        elemDG = reshape(1:3*NT,NT,3);
        nodeDG = node(elem(:),:);        
        elem2edge = Th.elem2edge;
        uhDG = zeros(3*NT,1);
        uhDG(1:NT) = -uh(elem2edge(:,1)) + uh(elem2edge(:,2)) + uh(elem2edge(:,3));
        uhDG(NT+(1:NT)) = uh(elem2edge(:,1)) - uh(elem2edge(:,2)) + uh(elem2edge(:,3));
        uhDG(2*NT+(1:NT)) = uh(elem2edge(:,1)) + uh(elem2edge(:,2)) - uh(elem2edge(:,3));
        showresult(nodeDG,elemDG,pde.exactu,uhDG);
        pause(1);
    end
    % compute error
    ErrL2(k) = getL2error(node,elem,pde.exactu,uh);
    ErrH1(k) = getH1error(node,elem,pde.Du,uh);
end

%% Plot convergence rates and display error table
figure(2);
showrateh(h,ErrH1,ErrL2);
fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e');

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (1st order) and L2-norm
% (2nd order) is observed. 
