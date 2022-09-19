clc;clear;close all

%% Parameters
maxIt = 5;
N = zeros(maxIt,1);    
h = zeros(maxIt,1);
Err = zeros(maxIt,1); 

%% Generate an intitial mesh
[node,elem] = squaremesh([0 1 0 1],1/5);
bdStr = [];

%% Exact eigenvalues
nm = [1 1; 
      1 2; 2 1;
      1 3; 2 2; 3 1;
      1 4; 2 3; 3 2; 4 1];
lamExact = pi^2*(nm(:,1).^2 + nm(:,2).^2);  
lamExact = sort(lamExact);  % in ascending order
nlam = length(lamExact);

%% Finite element method
Vh = 'P3'; quadOrder = 7;
for k = 1:maxIt
    % refine mesh
    [node,elem] = uniformrefine(node,elem);
    % get the mesh information
    Th = FeMesh2d(node,elem,bdStr);
    % A
    Coef  = 1;  Test  = 'v.grad';  Trial = 'u.grad';
    A = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
    NNdof = size(A,1);
    % B
    Coef  = 1;  Test  = 'v.val';  Trial = 'u.val';
    B = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
    % apply Dirichlet boundary conditions
    gDLogic = 1; % vector case: [u1,u2,u3] = [1,0,1]
    on = 1;
    [A,freedof] = apply2dMat(A,on,Th,Vh,gDLogic);
    B = apply2dMat(B,on,Th,Vh,gDLogic);

    % find eigenvalues
    Uh = zeros(NNdof,nlam);
    [Uh(freedof,:),D] = eigs(A,B,nlam,'smallestabs');
    lam = diag(D);

    % record
    N(k) = size(elem,1);
    h(k) = 1/(sqrt(size(node,1))-1);

    % compute errors
    Err(k) = norm(lam-lamExact)/sqrt(nlam);
end

%% Plot convergence rates and display results
figure,
showrateh(h,Err,'|\lambda - \lambda_h|');

fprintf('\n');
colname = {'Exact eigenvalues',  'Numerical eigenvalues'};
disptable(colname, lamExact,'%0.5f', lam,'%0.5f');