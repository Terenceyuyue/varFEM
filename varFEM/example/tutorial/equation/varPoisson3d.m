function uh = varPoisson3d(Th,pde,Vh,quadOrder)
%varPoisson3 Poisson equation in 3D: Pk Lagrange element (k<=2)

%   This function produces the finite element approximation of the 
%   Poisson equation
% 
%       -div(a*grad(u)) + cu = f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Robin boundary condition     g_R*u + a*grad(u)*n=g_N on \Gamma _R
%

% Quadrature orders for int1d, int2d and int3d
if nargin==2, Vh = 'P1'; quadOrder = 3; end % default: P1
if nargin==3, quadOrder = 3; end

%% Assemble stiffness matrix 
% Omega
Coef  = {pde.a, pde.c};
Test  = {'v.grad', 'v.val'};
Trial = {'u.grad', 'u.val'};
kk = assem3d(Th,Coef,Test,Trial,Vh,quadOrder);

% Gamma_R
if ~isempty(Th.bdStr)
    Th.on = 1;
    Coef  = pde.g_R;
    Test  = 'v.val';
    Trial = 'u.val'; 
    kk = kk + assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
end

%% Assemble the right hand side 
% Omega
Coef = pde.f;  Test = 'v.val';
ff = assem3d(Th,Coef,Test,[],Vh,quadOrder);
% Gamma_R
if ~isempty(Th.bdStr)
    %Coef = @(p) pde.g_R(p).*pde.uexact(p) + pde.a(p).*(pde.Du(p)*n');

    fun = @(p) pde.g_R(p).*pde.uexact(p);
    Cmat1 = interpFaceMat(fun,Th,quadOrder);
    fun = @(p) repmat(pde.a(p),1,3).*pde.Du(p);
    Cmat2 = interpFaceMat(fun,Th,quadOrder);
    Coef = Cmat1 + Cmat2;

    ff = ff + assem2d(Th,Coef,Test,[],Vh,quadOrder);
end

%% Apply Dirichlet boundary value conditions
g_D = pde.g_D;
on = 2 - 1*isempty(Th.bdStr);
uh = apply3d(on,Th,kk,ff,Vh,g_D);