clc;clear;close all
% ----------- Mesh and boudary conditions ----------
g = [2  2  2  2  2  2   % decomposed geometry matrix
     0  1  1 -1 -1  0
     1  1 -1 -1  0  0
     0  0  1  1 -1 -1
     0  1  1 -1 -1  0
     1  1  1  1  1  1
     0  0  0  0  0  0];
[p,e,t] = initmesh(g,'hmax',1); % initial mesh

bdNeumann = []; % only Dirichlet condition for elasticity3

% ------------------------ PDE data ------------------------
lambda = 1; mu = 1;
para.lambda = lambda; para.mu = mu;
pde = elasticitydata(para);

% ----------------- elasticity1 ---------------------
maxIt = 5;
N = zeros(maxIt,1);  h = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);  ErrH1 = zeros(maxIt,1);
errwL2 = zeros(maxIt,1);  errwH1 = zeros(maxIt,1);
for k = 1:maxIt
    [p,e,t] = refinemesh(g,p,e,t); % refine mesh
    node = p'; elem = t(1:3,:)';
    bdStruct = setboundary(node,elem,bdNeumann);
    Th.node = node; Th.elem = elem; Th.bdStruct = bdStruct;
    uh = elasticity3_variationalBlock(Th,pde);
    uh = reshape(uh,[],2);
    NT = size(elem,1);    h(k) = 1/sqrt(NT);
    
    tru = eye(2); trDu = eye(4);
    errL2 = zeros(1,2);  errH1 = zeros(1,2); % square
    for id = 1:2
        uid = uh(:,id);
        u = @(pz) pde.uexact(pz)*tru(:, id);
        Du = @(pz) pde.Du(pz)*trDu(:, 2*id-1:2*id);
        errL2(:,id) = getL2error(node,elem,uid,u);
        errH1(:,id) = getH1error(node,elem,uid,Du);
    end
    
    ErrL2(k) = sqrt(sum(errL2.^2,2));
    ErrH1(k) = sqrt(sum(errH1.^2,2));
end

% ---------- Plot convergence rates -----------
figure;
showrateh(h, ErrL2, ErrH1);
