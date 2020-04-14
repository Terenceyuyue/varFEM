function u = FEM1D(node,elem1D,pde,bdStruct)
%FEM1D solves the 1-D partial diffential equation 
%  
%  -au''+bu+cu = f,  x \in (x0,xL);
%     g_N = du(x0), g_D = u(xL);
%  or g_D = u(x0),  g_N = du(xL);
%
% where, g_N for Neumann boundary conditions and 
%        g_D for Dirichlet boundary conditions,
%        a, b and c are constants.
%
% Copyright (C) Terence YUE Yu. 

N = size(node,1); nel = size(elem1D,1); Ndof = 2;
% -------- Sparse assembling indices -----------
nnz = nel*Ndof^2;
ii = zeros(nnz,1); jj = zeros(nnz,1); 
id = 0; 
for i = 1:Ndof
    for j = 1:Ndof
        ii(id+1:id+nel) = elem1D(:,i);   % zi
        jj(id+1:id+nel) = elem1D(:,j);   % zj
        id = id + nel; 
    end
end

% ----------- Stiffness matrix -----------
% All element matrices 
para = pde.para;
a = para.a; b = para.b; c = para.c;
h = diff(node);
k11 = a./h+b/2*(-1)+c*h./6*2;      
k12 = a./h*(-1)+b/2+c*h./6;
k21 = a./h*(-1)+b/2*(-1)+c*h./6;   
k22 = a./h+b/2+c*h./6*2;
K = [k11,k12,k21,k22]; % stored in rows
% stiffness matrix
kk = sparse(ii,jj,K(:),N,N);

% ------------- Load vector ------------
x1 = node(1:N-1); x2 = node(2:N);
xc = (x1+x2)./2; 
F1 = pde.f(xc).*h./2; F2 = F1; F = [F1,F2];
ff = accumarray(elem1D(:), F(:),[N 1]);

% ------- Neumann boundary conditions -----------
Neumann = bdStruct.Neumann; 
if ~isempty(Neumann)
    nvec = 1;
    if find(elem1D(:,1)==Neumann), nvec = -1; end
    Dnu = pde.Du(node(Neumann,:))*nvec;
    ff(Neumann) = ff(Neumann) + a*Dnu;
end

% ------- Dirichlet boundary conditions ----------
Dirichlet = bdStruct.Dirichlet; g_D = pde.g_D;
isBdNode = false(N,1); isBdNode(Dirichlet) = true;
bdNode = (isBdNode); freeNode = (~isBdNode);
u = zeros(N,1); u(bdNode) = g_D(node(Dirichlet));
ff = ff - kk*u;

% ------------------ Solver -------------------
u(freeNode) = kk(freeNode,freeNode)\ff(freeNode);

