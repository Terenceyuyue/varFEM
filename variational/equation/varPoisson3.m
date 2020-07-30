function uh = varPoisson3(Th,pde,Vh,quadOrder)
%% VARPOISSON3 Poisson equation in 3D: Pk Lagrange element (k<=2)

%   This function produces the finite element approximation of the 
%   Poisson equation
% 
%       -div(a*grad(u)) + cu = f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Robin boundary condition     g_R*u + a*grad(u)*n=g_N on \Gamma _R

% Quadrature orders for int1d, int2d and int3d
if nargin==2, Vh = 'P1'; quadOrder = 3; end % default: P1
if nargin==3, quadOrder = 3; end

%% Get 2D mesh in boundary integrals
% elem2D = Gamma_R
bdStruct = Th.bdStruct;
Th.elem2D = bdStruct.bdFaceN; 

%% Assemble stiffness matrix 
% Omega
Coef  = {pde.a, pde.c};
Test  = {'v.grad', 'v.val'};
Trial = {'u.grad', 'u.val'};
kk = int3d(Th,Coef,Test,Trial,Vh,quadOrder);

% Gamma_R
if ~isempty(Th.elem2D)
    Coef  = {pde.g_R};
    Test  = {'v.val'};
    Trial = {'u.val'}; 
    kk = kk + int2d(Th,Coef,Test,Trial,Vh,quadOrder);
end

%% Assemble the right hand side 
% Omega
Coef = pde.f;  Test = 'v.val';
ff = int3d(Th,Coef,Test,[],Vh,quadOrder);
% Gamma_R
if ~isempty(Th.elem2D)
    %Coef = @(p) g_R(p).*pde.uexact(p) + pde.a(p).*(pde.Du(p)*n');
    Cmat_gR = getMat2d(pde.g_R,Th,quadOrder);
    Cmat_u = getMat2d(pde.exactu,Th,quadOrder);
    Cmat_a = getMat2d(pde.a,Th,quadOrder);
    Cmat_Dnu = getMat2d(pde.Du,Th,quadOrder);
    Coef = Cmat_gR.*Cmat_u + Cmat_a.*Cmat_Dnu;
    ff = ff + int2d(Th,Coef,Test,[],Vh,quadOrder);
end


%% Apply Dirichlet boundary value conditions
g_D = pde.g_D;
uh = Applyboundary3D(Th,kk,ff,g_D,Vh);