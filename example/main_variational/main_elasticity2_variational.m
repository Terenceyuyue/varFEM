clc;clear;close all
tic;
% --------------- Mesh and boudary conditions ---------------
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
Nx = 4; Ny = 4; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);

bdNeumann = 'abs(y-0)<1e-4 | abs(x-1)<1e-4'; % string for Neumann

% ------------------------ PDE data ------------------------
lambda = 1; mu = 1;
para.lambda = lambda; para.mu = mu;
pde = elasticitydata(para);

% ----------------- elasticity1 ---------------------
maxIt = 5;
N = zeros(maxIt,1);  h = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);  ErrH1 = zeros(maxIt,1);
errwL2 = zeros(maxIt,1);  errwH1 = zeros(maxIt,1);

feSpace = 'P1';
if strcmpi(feSpace,'P1'), quadOrder = 3; end
if strcmpi(feSpace,'P2'), quadOrder = 4; end
if strcmpi(feSpace,'P3'), quadOrder = 5; end
for k = 1:maxIt
    [node,elem] = uniformrefine(node,elem);
    bdStruct = setboundary(node,elem,bdNeumann);
    Th = getTh(node,elem); 
    uh = elasticity2_variational(Th,pde,feSpace,quadOrder);
    uh = reshape(uh,[],2);
    NT = size(elem,1);    h(k) = 1/sqrt(NT);
    
    tru = eye(2); trDu = eye(4);
    errL2 = zeros(1,2);  errH1 = zeros(1,2); % square
    for id = 1:2
        uid = uh(:,id);
        u = @(pz) pde.uexact(pz)*tru(:, id);
        Du = @(pz) pde.Du(pz)*trDu(:, 2*id-1:2*id);
        errL2(:,id) = getL2error(node,elem,uid,u,feSpace,quadOrder);
        errH1(:,id) = getH1error(node,elem,uid,Du,feSpace,quadOrder);
    end
    
    ErrL2(k) = sqrt(sum(errL2.^2,2));
    ErrH1(k) = sqrt(sum(errH1.^2,2));
end

% ---------- Plot convergence rates -----------
figure;
showrateh(h, ErrL2, ErrH1);

hold on
toc