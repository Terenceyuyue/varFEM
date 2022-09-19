function cc = interp2dMat(uh,intType,Th,Vh,quadOrder)
% uh is the solution vector of FEM or a function handle
% cc is the matrix for numerical integration (see Coef matrix in
% assem2d.m)
% intType: .val, .dx, .dy or .grad, ...

node = Th.node; elem = Th.elem;
NT = size(elem,1);
if isa(uh, 'function_handle')
    if mycontains(intType,'.val')
        [lambda,weight] = quadpts(quadOrder);
        ng = length(weight);
        cc = zeros(NT,ng);
        z1 = node(elem(:,1),:);
        z2 = node(elem(:,2),:);
        z3 = node(elem(:,3),:);
        for p = 1:ng
            pz = lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
            cc(:,p) = uh(pz);
        end
        return;
    end
    uh = interp2d(uh,Th,Vh); % for other cases, use FEM approximation
end

%% Basis function
phi = Base2d(intType,node,elem,Vh,quadOrder); % phi1,phi2,phi3
elem2dof = dof2d(Th,Vh);

%% Integration matrix
nG = size(phi{1},2); % nG = 2*ng for .grad
cc = zeros(NT,nG);
for ib = 1:length(phi)
    ui = uh(elem2dof(:,ib));
    cc = cc + repmat(ui,1,nG).*phi{ib};   % ui*phi_i
end