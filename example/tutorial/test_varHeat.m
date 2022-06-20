%%test_varHeat: Pk Lagrange element (k<=3)
%
%   Heat equation:
%
%   D.E.  u_t - \Delta u = f  in \Omega
%   I.C.  u(x,y,0) = u0(x,y)  
%   B.C.  u = g_D  on \Gamma_D       (Dirichlet)
%         grad(u)*n=g_N on \Gamma_N  (Neumann)
%

clc;clear;close all

%% Mesh
% generate mesh
Nx = 10;
[node,elem] = squaremesh([0 1 0 1],1/Nx);
% mesh info
bdStr = 'x==0'; % Neumann
Th = FeMesh2d(node,elem,bdStr);
% time
Nt = Nx^2;
t = linspace(0,1,Nt+1)';  dt = t(2)-t(1);

%% PDE data
pde = heatData();

%% Bilinear form
Vh = 'P1';  quadOrder = 7;
Coef  = {1/dt, 1};
Test  = {'v.val', 'v.grad'};
Trial = {'u.val', 'u.grad'};
kk = assem2d(Th,Coef,Test,Trial,Vh,quadOrder); 

%% Backward Euler
u0 = @(p) pde.uexact(p,t(1));
uh0 = interp2d(u0,Th,Vh); % dof vector
uf = zeros(Nt+1,2);  % record solutions at p-th point
p = 2*Nx;
uf(1,:) = [uh0(p),uh0(p)]; 
for n = 1:Nt
    % Linear form
    fun = @(p) pde.f(p, t(n+1));
    Coef = fun;
    Test = 'v.val';    
    ff = assem2d(Th,Coef,Test,[],Vh,quadOrder);
    Coef = uh0/dt;
    ff = ff + assem2d(Th,Coef,Test,[],Vh,quadOrder);

    % Neumann boundary condition
    if ~isempty(bdStr)
        Th.on = 1;
        fun = @(p) pde.Du(p,t(n+1));
        Coef = interpEdgeMat(fun,Th,quadOrder);
        %[fh,Nh] = interpEdge(fun,Th,Vh);
        %Coef = sum(fh.*Nh,2);

        ff = ff + assem1d(Th,Coef,Test,[],Vh,quadOrder);
    end

    % Dirichlet boundary condition
    ue = @(p) pde.uexact(p,t(n+1));
    on = 2 - 1*isempty(bdStr); % on = 1 for the whole boundary if bdStr = []
    uh = apply2d(on,Th,kk,ff,Vh,ue);

    % Error
    if rem(n,10)==0
        ErrL2 = varGetL2Error(Th,ue,uh,Vh,quadOrder);
        fprintf('t = %.4f | ErrL2 = %.4e \n', t(n+1), ErrL2);
    end

    % Record
    uhe = interp2d(ue,Th,Vh);
    uf(n+1,:) = [uhe(p), uh(p)];

    % Update
    uh0 = uh;
end

figure,
showresult(node,elem,ue,uh);

figure,
plot(t,uf(:,1),'r-*',t,uf(:,2),'k-','linewidth',1);
xlabel('t');
