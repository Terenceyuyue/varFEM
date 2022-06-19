function [uh,ph] = varNavierStokes(Th,pde,Vh,quadOrder)
%varNavierStokes Taylor-Hood (P2-P1) FEM for steady Navier-Stokes problem
%       u = [u1, u2]
%       -nu*\Delta u + (u.grad)u + grad p = f   in \Omega,
%       -div u = 0  in  \Omega,
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma
%

% Quadrature orders for int1d and int2d
if nargin==2
    Vh = {'P2','P2','P1'}; quadOrder = 5;
end % default: Taylor-Hood
if nargin==3, quadOrder = 5; end

nu = pde.nu;  eps = 1e-10;
vstr = {'v1','v2','q'}; ustr = {'u1','u2','p'};

%% Fixed bilinear forms
% a(uh^{n+1), v) +  b(v, ph^{n+1}) + b(uh^{n+1},q) - eps*p*q
Coef = { nu, nu,  -1, -1,  -1, -1, -eps};
Test  = {'v1.grad','v2.grad',  'v1.dx','v2.dy',  'q.val', 'q.val', 'q.val'};
Trial = {'u1.grad','u2.grad',  'p.val','p.val',  'u1.dx', 'u2.dy', 'p.val'};

[Test,Trial] = getStdvarForm(vstr, Test, ustr, Trial);
kk = int2d(Th,Coef,Test,Trial,Vh,quadOrder);

%% Fixed linear forms
f1 = @(p) pde.f(p)*[1;0];
f2 = @(p) pde.f(p)*[0;1];
Coef = {f1, f2};
Test = {'v1.val', 'v2.val'};
ff = int2d(Th,Coef,Test,[],Vh,quadOrder);

%% Newton iteration
% initial data
u1 = @(p) 0*p(:,1);
u2 = @(p) 0*p(:,1);
p = @(p) 0*p(:,1);
uh1 = interp2d(u1,Th,Vh{1});
uh2 = interp2d(u2,Th,Vh{2});
ph = interp2d(p,Th,Vh{3});

erru = 1;  n = 1; 
maxIt = 20;
while n<maxIt && erru>1e-10

    % transform FEM functions to Coef matrices
    u1c = interp2dMat(uh1,'u1.val',Th,Vh{1},quadOrder);
    u2c = interp2dMat(uh2,'u2.val',Th,Vh{2},quadOrder);

    u1xc = interp2dMat(uh1,'u1.dx',Th,Vh{1},quadOrder);
    u1yc = interp2dMat(uh1,'u1.dy',Th,Vh{1},quadOrder);
    u2xc = interp2dMat(uh2,'u2.dx',Th,Vh{2},quadOrder);
    u2yc = interp2dMat(uh2,'u2.dy',Th,Vh{2},quadOrder);

    % bilinear form: N(uh^{n+1}; uh^n, v) + N(uh^n; uh^{n+1}, v)
    Coef = {u1xc, u1yc, u2xc, u2yc, ...
        u1c, u2c, u1c, u2c};
    Test = {'v1.val','v1.val','v2.val','v2.val', ...
        'v1.val','v1.val','v2.val','v2.val'};
    Trial = {'u1.val','u2.val','u1.val','u2.val', ...
        'v1.dx','v1.dy','v2.dx','v2.dy'};
    [kkN,info] = int2d(Th,Coef,Test,Trial,Vh,quadOrder);
    kkn = kk + kkN;

    % linear form: N(uh^n; uh^n, v)
    Coef = {u1c.*u1xc+u2c.*u1yc,  u1c.*u2xc+u2c.*u2yc};
    Test = {'v1.val', 'v2.val'};
    ffn = ff + int2d(Th,Coef,Test,[],Vh,quadOrder);

    % Dirichlet boundary conditions
    g_D1 = @(p) pde.g_D(p)*[1;0];
    g_D2 = @(p) pde.g_D(p)*[0;1];
    g_D = {g_D1, g_D2, []};
    on = 1;
    U = apply2d(on,Th,kkn,ffn,Vh,g_D);
    
    % Update
    NNdofu = info.NNdofu;
    id1 =  NNdofu(1);  id2 = NNdofu(1)+ NNdofu(2);
    erru = norm([U(1:id1)-uh1; U(id1+1:id2)-uh2]);

    uh1 = U(1:id1);
    uh2 = U(id1+1:id2);
    ph = U(id2+1:end);
    n = n+1;
end
uh = [uh1,uh2];

fprintf('Number of iterations = %d, \t  erru = %.4e \n', n,erru);