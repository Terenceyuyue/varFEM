function u = Applyboundary1D(Th,kk,ff,pde)

% NNdof
node = Th.node;  NNdof = size(ff,1);

% boundary information
Coef = pde.a;  
bdStruct = Th.bdStruct;  
Dirichlet = bdStruct.Dirichlet; Neumann = bdStruct.Neumann;

% ------- Neumann boundary conditions -----------
if ~isempty(Neumann)
    nvec = 1;
    if find(Th.elem1D(:,1)==Neumann), nvec = -1; end
    Dnu = pde.Du(node(Neumann,:))*nvec;
    ff(Neumann) = ff(Neumann) + Coef(node(Neumann,:))*Dnu;
end

% ------- Dirichlet boundary conditions ----------
isBdNode = false(NNdof,1); isBdNode(Dirichlet) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
u = zeros(NNdof,1); u(bdDof) = pde.g_D(node(Dirichlet,:));
ff = ff - kk*u;

% ------------------ Solver ----------------
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);