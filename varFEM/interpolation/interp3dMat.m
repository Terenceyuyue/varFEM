function cc = interp3dMat(uh,intType,Th,Vh,quadOrder)
% uh is the solution vector of FEM or a function handle
% cc is the matrix for numerical integration (see Coef matrix in
% assem3d.m)
% intType: .val, .dx, .dy, dz or .grad

node = Th.node; elem3 = Th.elem3;  
NT = size(elem3,1);
if isa(uh, 'function_handle')
    [lambda,weight] = quadpts3(quadOrder);
    ng = length(weight);  
    cc = zeros(NT,ng);
    z1 = node(elem3(:,1),:);
    z2 = node(elem3(:,2),:);
    z3 = node(elem3(:,3),:);
    z4 = node(elem3(:,4),:);
    for p = 1:ng
        pz = lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3 + lambda(p,4)*z4;
        cc(:,p) = uh(pz);
    end
    return;
end

%% Basis function
phi = Base3d(intType,node,elem3,Vh,quadOrder); % phi1,phi2,phi3
elem2dof = dof3d(Th,Vh);

%% Integration matrix
nG = size(phi{1},3); % nG = 3*ng for .grad
cc = zeros(NT,nG);
for ib = 1:length(phi)
    ui = uh(elem2dof(:,ib));
    cc = cc + repmat(ui,1,nG).*phi{ib};   % ui*phi_i
end