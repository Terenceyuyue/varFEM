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
for i = 1:3
    [p,e,t] = refinemesh(g,p,e,t); % refine mesh
end
node = p'; elem = t(1:3,:)';
figure,showmesh(node,elem);
bdStruct = setboundary(node,elem);

% ------------ PDE data ------------
lambda = 1; mu = 1;
para.lambda = lambda; para.mu = mu;
pde = elasticitydata1(para);

% ----------- elasticity --------
u = elasticityLaplace(node,elem,pde,bdStruct);

% --------- error analysis -------
u = reshape(u,[],2);
uexact = pde.uexact;  ue = uexact(node);
id = 1;
figure,
subplot(1,2,1), showsolution(node,elem,u(:,id));
zlabel('u');
subplot(1,2,2), showsolution(node,elem,ue(:,id));
zlabel('ue');
Eabs = u-ue;  % Absolute errors
figure,showsolution(node,elem,Eabs(:,id)); zlim('auto');
format shorte
Err = norm(Eabs)./norm(ue)
