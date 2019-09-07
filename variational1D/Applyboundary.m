function u = Applyboundary(Th,kk,ff,acoef,pde,bdStruct)

node = Th.node; N = length(node); 
if isa(acoef,'double'), acoef = @(x) acoef+0*x; end

% ------- Neumann boundary conditions -----------
Neumann = bdStruct.Neumann;
if ~isempty(Neumann)
    bn = [-1;zeros(N-2,1);1]; % -1: left norm vector
    bn(Neumann) = bn(Neumann)*pde.g_N(node(Neumann));
    ff = ff - acoef(node(Neumann))*bn;
end

% ------- Dirichlet boundary conditions ----------
Dirichlet = bdStruct.Dirichlet;
isBdNode = false(N,1); isBdNode(Dirichlet) = true;
bdNode = (isBdNode); freeNode = (~isBdNode);
u = zeros(N,1); u(bdNode) = pde.g_D(node(bdNode));
ff = ff - kk*u;

% ------------------ Solver ----------------
u(freeNode) = kk(freeNode,freeNode)\ff(freeNode);