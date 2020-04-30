function u = Poisson_variational(Th,pde,Vh,quadOrder)

% Quadrature orders for int1d and int2d
if nargin==2, Vh = 'P1'; quadOrder = 3; end % default: P1
if nargin==3, quadOrder = 3; end

% ------------------------ Mesh Th --------------------------
% elem1D associated with Gamma_R
bdStruct = Th.bdStruct;
Th.elem1D = bdStruct.elemN; Th.bdIndex1D = bdStruct.bdIndexN;

% ---------------------- Stiffness matrix -----------------------
% Omega
Coef  = {pde.a, pde.c};
Test  = {'v.grad', 'v.val'};
Trial = {'u.grad', 'u.val'};
kk = int2d(Th,Coef,Test,Trial,Vh,quadOrder);

% Gamma_R
if ~isempty(Th.elem1D)
    Coef  = {pde.g_R};
    Test  = {'v.val'};
    Trial = {'u.val'};    
    kk = kk + int1d(Th,Coef,Test,Trial,Vh,quadOrder);
end

% -------------------------- Load vector -------------------------
% Omega
Coef = pde.f;  Test = 'v.val';
ff = int2d(Th,Coef,Test,[],Vh,quadOrder);
% Gamma_R
if ~isempty(Th.elem1D)
    %Coef = @(p) g_R(p).*pde.uexact(p) + pde.a(p).*(pde.Du(p)*n');
    Cmat_gR = getMat1d(pde.g_R,Th,quadOrder);
    Cmat_u = getMat1d(pde.uexact,Th,quadOrder);
    Cmat_a = getMat1d(pde.a,Th,quadOrder);
    Cmat_Dnu = getMat1d(pde.Du,Th,quadOrder);
    Coef = Cmat_gR.*Cmat_u + Cmat_a.*Cmat_Dnu;
    ff = ff + int1d(Th,Coef,Test,[],Vh,quadOrder);
end

% ------------------- Boundary value conditions ------------------
g_D = pde.g_D;
u = Applyboundary2D(Th,kk,ff,g_D,Vh);
