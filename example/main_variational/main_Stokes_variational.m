clc;clear;close all
tic;
% --------------- Mesh and boudary conditions ---------------
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
Nx = 2; Ny = 2; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);

% ------------------------ PDE data ------------------------
pde = Stokesdata();

% ----------------- elasticity1 ---------------------
maxIt = 5;
N = zeros(maxIt,1);  h = zeros(maxIt,1);
ErruL2 = zeros(maxIt,1);  ErruH1 = zeros(maxIt,1);
ErrpL2 = zeros(maxIt,1);  ErrpH1 = zeros(maxIt,1);

Vh = {'P2','P2','P1'}; % v = [v1,v2,q] --> [v1,v2,v3]
quadOrder = 4;
for k = 1:maxIt
    [node,elem] = uniformrefine(node,elem);
    Th = getTh(node,elem); 
    [uh,ph] = Stokes_variational(Th,pde,Vh,quadOrder);
    NT = size(elem,1);    h(k) = 1/sqrt(NT);
    
    tru = eye(2); trDu = eye(4);
    erruL2 = zeros(1,2);  erruH1 = zeros(1,2); % square
    for id = 1:2
        uid = uh(:,id);
        u = @(pz) pde.uexact(pz)*tru(:, id);
        Du = @(pz) pde.Du(pz)*trDu(:, 2*id-1:2*id);
        erruL2(:,id) = getL2error(node,elem,uid,u,Vh{1},4);
        erruH1(:,id) = getH1error(node,elem,uid,Du,Vh{1},4);
    end
    
    ErruL2(k) = sqrt(sum(erruL2.^2,2));
    ErruH1(k) = sqrt(sum(erruH1.^2,2));
    
    ErrpL2(k) = getL2error(node,elem,ph,pde.pexact,Vh{3},4);
    ErrpH1(k) = getH1error(node,elem,ph,pde.Dpexact,Vh{3},4);
end

% ---------- Plot convergence rates -----------
figure;
showrateh(h, ErruL2, ErruH1, '||u - u_h||', '||Du - Du_h||');
figure;
showrateh(h, ErrpL2, ErrpH1, '||p - p_h||', '||Dp - Dp_h||');

toc