function u = apply1d(Th,kk,ff,pde)
%%apply1d deals with Neumann and Dirichlet boundary conditions of 1-D problems
%

% Mesh
node = Th.node;  elem1d = Th.elem1d;
NNdof = length(ff);

%% Neumann boundary condition
bdStr = Th.bdStr;
if ~isempty(bdStr)
    Neumann = Th.bdNodeIdxType{1};
    nvec = 1;
    for i = 1:length(Neumann)
        if find(elem1d(:,1)==Neumann(i)), nvec = -1; end
        Dnu = pde.Du(node(Neumann(i),:))*nvec;
        ff(Neumann(i)) = ff(Neumann(i)) + pde.a(node(Neumann(i),:))*Dnu;
    end
end

%% Dirichlet boundary condition
if isempty(bdStr)  % a and b are Dirichlet
    on = 1; 
end
if length(Th.bdNodeIdxType)==2 % a or b is Dirichlet
    on = 2; 
end
if ~isempty(bdStr) && length(Th.bdNodeIdxType)==1 % no Dirichlet
    %on = [];
    error('Dirichlet condition is necessary');
end
Dirichlet = Th.bdNodeIdxType{on};
isBdDof = false(NNdof,1); isBdDof(Dirichlet) = true;
bddof = (isBdDof);  freedof = (~isBdDof);
bdval = pde.g_D(node(Dirichlet,:));
u = zeros(NNdof,1);  u(bddof) = bdval;
ff = ff - kk*u;

%% Solver
u(freedof) = kk(freedof,freedof)\ff(freedof);