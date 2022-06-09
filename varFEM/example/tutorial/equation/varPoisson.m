function uh = varPoisson(Th,pde,Vh,quadOrder)
%varPoisson Poisson equation: Pk Lagrange element (k<=3)

%   This function produces the finite element approximation of the 
%   Poisson equation
% 
%       -div(a*grad(u)) + cu = f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Robin boundary condition     g_R*u + a*grad(u)*n=g_N on \Gamma _R
%

% Quadrature orders for assem1d and assem2d
if nargin==2, Vh = 'P1'; quadOrder = 3; end % default: P1
if nargin==3, quadOrder = 3; end

%% Assemble stiffness matrix 
% Omega
Coef  = {pde.a, pde.c};
Test  = {'v.grad', 'v.val'};
Trial = Test; %{'u.grad', 'u.val'};
kk = int2d(Th,Coef,Test,Trial,Vh,quadOrder); % or assem2d

% Robin boundary condition
bdStr = Th.bdStr;
if ~isempty(bdStr)
    Th.elem1d = Th.bdEdgeType{1};
    Th.elem1dIdx = Th.bdEdgeIdxType{1};

    Coef  = {pde.g_R};
    Test  = {'v.val'};
    Trial = {'u.val'};
    kk = kk + int1d(Th,Coef,Test,Trial,Vh,quadOrder); % or assem1d
end

%% Assemble the right hand side 
% Omega
Coef = pde.f;
Test = 'v.val';
ff = int2d(Th,Coef,Test,[],Vh,quadOrder);
% Neumann boundary condition
if ~isempty(bdStr)
    %Coef = @(p) pde.g_R(p).*pde.uexact(p) + pde.a(p).*(pde.Du(p)*n');

    fun = @(p) pde.g_R(p).*pde.uexact(p);
    Cmat1 = interpEdgeMat(fun,Th,quadOrder);
    fun = @(p) repmat(pde.a(p),1,2).*pde.Du(p);
    Cmat2 = interpEdgeMat(fun,Th,quadOrder);
    Coef = Cmat1 + Cmat2;

    ff = ff + int1d(Th,Coef,Test,[],Vh,quadOrder);
end

%% Apply Dirichlet boundary value conditions
g_D = pde.g_D;
on = 2 - 1*isempty(bdStr);
uh = apply2d(on,Th,kk,ff,g_D,Vh);