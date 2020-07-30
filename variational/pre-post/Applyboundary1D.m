function u = Applyboundary1D(Th,A,b,pde)
%% APPLYBOUNDARY1D deals with the Neumann and Dirichlet boundary conditions of 1D problems

% NNdof
node = Th.node;  elem = Th.elem;
NNdof = length(b);

% boundary information
Coef = pde.a;  
bdStruct = Th.bdStruct;
Dirichlet = bdStruct.Dirichlet; 
Neumann = bdStruct.Neumann;

% Neumann boundary conditions
if ~isempty(Neumann)
    nvec = 1;
    if find(elem(:,1)==Neumann), nvec = -1; end
    Dnu = pde.Du(node(Neumann,:))*nvec;
    b(Neumann) = b(Neumann) + Coef(node(Neumann,:))*Dnu;
end

% Dirichlet boundary conditions
isBdDof = false(NNdof,1); isBdDof(Dirichlet) = true;
bddof = (isBdDof);  freedof = (~isBdDof);
bdval = pde.g_D(node(Dirichlet,:));
u = zeros(NNdof,1);  u(bddof) = bdval;
b = b - A*u;

%% Solver
solver = 'amg';
if NNdof < 2e3, solver = 'direct'; end
% solve
switch solver
    case 'direct'
        u(freedof) = A(freedof,freedof)\b(freedof);
    case 'amg'
        option.solver = 'CG';
        u(freedof) = amg(A(freedof,freedof),b(freedof),option);                 
end

