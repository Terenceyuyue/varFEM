clc;clear;close all
%%% This is an exmple given in FreeFem Documentation: Release 4.6
%   Subsection 2.12 - Newton Method for the Steady Navier-Stokes equations 
%   (see NSNewton.edp)
%
% We check the main step of defining the linear system in this M-file.
% The postprocessing for the blowup of the Newton method will be treated in
% test_NSNewton.m
%

Vh = {'P2','P2','P1'}; % v = [v1,v2,q] --> [v1,v2,v3]
quadOrder = 7;

%% Mesh
[node,elem] = getMeshFreeFEM('meshdata_NSNewton.msh');
%load meshdata_NSNewton
%showmesh(node,elem);
% mesh info
bdStr = 'x==15'; % 1-Neumann
Th = FeMesh2d(node,elem,bdStr);

%% Assemble stiffness matrix
vstr = {'v1','v2','q'}; ustr = {'du1','du2','dp'};
eps = 1e-8;  nu = 1/50;

% initial data
u1 = @(p) double( (p(:,1).^2 + p(:,2).^2 > 2) );
u2 = @(p) 0*p(:,1);
p = @(p) 0*p(:,1);

uh1 = interp2d(u1,Th,Vh{1});
uh2 = interp2d(u2,Th,Vh{2});
ph = interp2d(p,Th,Vh{3});

for n = 1:5

    u1c = interp2dMat(uh1,'u1.val',Th,Vh{1},quadOrder);
    u2c = interp2dMat(uh2,'u2.val',Th,Vh{2},quadOrder);
    pc = interp2dMat(ph,'p.val',Th,Vh{3},quadOrder);

    u1xc = interp2dMat(uh1,'u1.dx',Th,Vh{1},quadOrder);
    u1yc = interp2dMat(uh1,'u1.dy',Th,Vh{1},quadOrder);
    u2xc = interp2dMat(uh2,'u2.dx',Th,Vh{2},quadOrder);
    u2yc = interp2dMat(uh2,'u2.dy',Th,Vh{2},quadOrder);

    Coef = { u1xc, u1yc, u2xc, u2yc,  ... % term 1
        u1c, u2c, u1c, u2c,  ... % term 2
        nu, nu, nu, nu, ... % term 3
        -1, -1, ... % term 4
        -1, -1, ... % term 5
        -eps ... % stablization term
        };
    Test  = { 'v1.val', 'v1.val', 'v2.val', 'v2.val', ... % term 1
        'v1.val', 'v1.val', 'v2.val', 'v2.val', ... % term 2
        'v1.dx', 'v1.dy', 'v2.dx', 'v2.dy', ... % term 3
        'v1.dx', 'v2.dy', ... % term 4
        'q.val', 'q.val', ... % term 5
        'q.val' ... % stablization term
        };
    Trial = { 'du1.val', 'du2.val', 'du1.val', 'du2.val',... % term 1
        'du1.dx', 'du1.dy', 'du2.dx', 'du2.dy', ... % term 2 
        'du1.dx', 'du1.dy', 'du2.dx', 'du2.dy', ... % term 3
        'dp.val', 'dp.val', ... % term 4
        'du1.dx', 'du2.dy', ... % term 5
        'dp.val' ... % stablization term
        };
    [Test,Trial] = getStdvarForm(vstr, Test, ustr, Trial); 
    [kk,info] = int2d(Th,Coef,Test,Trial,Vh,quadOrder);

    %% Assemble right hand side
    Coef = {u1c.*u1xc + u2c.*u1yc,  u1c.*u2xc + u2c.*u2yc, ...
        nu*u1xc, nu*u1yc, nu*u2xc, nu*u2yc, ...
        -pc, -pc, ...
        -(u1xc + u2yc),...
        -eps*pc   % stablization term
        };
    Test = {'v1.val', 'v2.val', ...
        'v1.dx', 'v1.dy', 'v2.dx', 'v2.dy', ...
        'v1.dx', 'v2.dy', ...
        'q.val',...
        'q.val'  % stablization term
        };
    Test = getStdvarForm(vstr, Test);
    ff = int2d(Th,Coef,Test,[],Vh,quadOrder);

    %% Apply Dirichlet boundary conditions
    g_D1 = @(p) 0*p(:,1);
    g_D2 = @(p) 0*p(:,1);
    g_D = {g_D1, g_D2, []};
    on = 2;
    U = apply2d(on,Th,kk,ff,Vh,g_D);
    NNdofu = info.NNdofu;
    id1 =  NNdofu(1);  id2 = NNdofu(1)+ NNdofu(2);
    duh1 = U(1:id1);
    duh2 = U(id1+1:id2);
    dph = U(id2+1:end);

    %% Update
    uh1 = uh1 - duh1;
    uh2 = uh2 - duh2;
    ph = ph - dph;
    
    %% Plot
    figure(1),
    clf;
    % surf
    subplot(2,2,1),
    showsolution(node,elem,ph);
    subplot(2,2,3)
    pff = solFreeFEM(sprintf('sol_NSNewtonp%d.txt',n));    
    showsolution(node,elem,pff);
    % contour
    x = [-5:0.05:-1, -1:0.01:1, 1:0.1:15];
    y = [-5:0.05:-1, -1:0.01:1, 1:0.05:5];
    subplot(2,2,2),
    varcontourf(x,y,node,elem,ph(1:Th.N),20); 
    title('Numerical solution of varFEM');
    subplot(2,2,4),
    varcontourf(x,y,node,elem,pff(1:Th.N),20); 
    title('Numerical solution of FreeFEM++');
    drawnow;
end