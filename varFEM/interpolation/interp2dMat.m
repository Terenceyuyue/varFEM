function cc = interp2dMat(uh,intType,Th,Vh,quadOrder)
% uh is the solution vector of FEM or a function handle
% cc is the matrix for numerical integration (see coef matrix in
% assem2d.m)
% intType: .val, .dx, .dy or .grad

if isa(uh, 'function_handle')
    uh = interp2d(uh,Th,Vh);
end

%% Basis function
phi = Base2d(intType,Th.node,Th.elem,Vh,quadOrder); % phi1,phi2,phi3
elem2dof = dof2d(Th,Vh);

%% Integration matrix
nG = size(phi{1},2); % nG = 2*ng for .grad
cc = zeros(Th.NT,nG);
for ib = 1:length(phi)
    ui = uh(elem2dof(:,ib));
    cc = cc + repmat(ui,1,nG).*phi{ib};   % ui*phi_i
end