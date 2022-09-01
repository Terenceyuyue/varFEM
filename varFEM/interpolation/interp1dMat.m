function cc = interp1dMat(uh,intType,Th,Vh,quadOrder)
% uh is the solution vector of FEM or a function handle
% cc is the matrix for numerical integration (see Coef matrix in
% assem1d.m)
% intType: .val, .dx

node = Th.node; elem1d = Th.elem1d;  
nel = size(elem1d,1);
if isa(uh, 'function_handle')
    [lambda,weight] = quadpts1(quadOrder);
    ng = length(weight);  
    cc = zeros(nel,ng);
    z1 = node(elem1d(:,1),:);
    z2 = node(elem1d(:,2),:);
    for p = 1:ng
        pz = lambda(p,1)*z1 + lambda(p,2)*z2;
        cc(:,p) = uh(pz);
    end
    return;
end

%% Basis function
phi = Base1d(intType,node,elem1d,Vh,quadOrder); 
elem2dof = dof1d(Th,Vh);

%% Integration matrix
nG = size(phi{1},2); 
cc = zeros(nel,nG);
for ib = 1:length(phi)
    ui = uh(elem2dof(:,ib));
    cc = cc + repmat(ui,1,nG).*phi{ib};   % ui*phi_i
end