function [uh,ErrH1,ErrL2] = varHeat(Th,pde,Vh,quadOrder)
%%varHeat: Pk Lagrange element (k<=3) with backward Euler
%
% Heat equation:
%
%   D.E.  u_t - \Delta u = f  in \Omega
%   I.C.  u(x,y,0) = u0(x,y)  
%   B.C.  u = g_D  on \Gamma_D       (Dirichlet)
%         grad(u)*n=g_N on \Gamma_N  (Neumann)
%
% The variational problem is to find u^n \in H_g^1(\Omega) such that
%
%    u^n - u^{n-1}
%  (---------------,  v)  +  (grad u^n, grad v) = (f^n, v) + (g_N, v)_{Gamma_N}  
%        dt 
%
% for n = 1,2,..., N_t-1,  where v \in H_0^1(\Omega)
%

% Quadrature orders for assem1d and assem2d
if nargin==2, Vh = 'P1'; quadOrder = 3; end % default: P1
if nargin==3, quadOrder = 3; end

%% Set time step
Nt = pde.Nt;
t = linspace(pde.t0, pde.tf, Nt+1)';  
dt = t(2)-t(1);

%% Assemble stiffness matrix 
Coef  = {1/dt, 1};
Test  = {'v.val', 'v.grad'};
Trial = Test; %{'u.val', 'u.grad'};
kk = assem2d(Th,Coef,Test,Trial,Vh,quadOrder); % stiffness matrix

%% Backward Euler
u0 = @(p) pde.uexact(p,t(1));
uh0 = interp2d(u0,Th,Vh);
bdStr = Th.bdStr;
n = 1;
while t(n)<pde.tn
    % Linear form
    fun = @(p) pde.f(p, t(n+1));
    Coef = fun;     Test = 'v.val';
    ff = assem2d(Th,Coef,Test,[],Vh,quadOrder);
    Coef = uh0/dt;  Test = 'v.val';  
    ff = ff + assem2d(Th,Coef,Test,[],Vh,quadOrder);

    % Neumann boundary condition
    Due = @(p) pde.Du(p,t(n+1));
    if ~isempty(bdStr)
        Th.on = 1;
        Coef = interpEdgeMat(Due,Th,quadOrder);
        %[fh,Nh] = interpEdge(Due,Th,Vh);
        %Coef = sum(fh.*Nh,2);
        ff = ff + assem1d(Th,Coef,Test,[],Vh,quadOrder);
    end

    % Dirichlet boundary condition
    ue = @(p) pde.uexact(p,t(n+1));
    on = 2 - 1*isempty(bdStr); % on = 1 for the whole boundary if bdStr = []
    uh = apply2d(on,Th,kk,ff,Vh,ue);    

    % Update
    uh0 = uh;
    n = n+1;
end

%% Compute error at the final step
ErrH1 = varGetH1Error(Th,Due,uh,Vh,quadOrder);
ErrL2 = varGetL2Error(Th,ue,uh,Vh,quadOrder);