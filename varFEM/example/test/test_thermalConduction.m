clc; clear; close all
%%% This is an exmple given in FreeFem++: Release 4.6
%   Subsection 2.6: Thermal Conduction (thermal.edp)

%% Mesh
% mesh generated by freefem++        
[node,elem] = getMeshFreeFEM('meshdata_thermalConduction.msh');

% mesh info
bdStr = 'y==0 | y==1';
Th = FeMesh2d(node,elem,bdStr);

%% Parameters
u0 = @(p) 10 + 90*p(:,1)/6;
k = @(p) 1.8*(p(:,2)<0.5) + 0.2;

ue = 25; % environment
alpha = 0.25; 
T = 5; dt = 0.1;

%% Bilinear form
Vh = 'P1';  quadOrder = 5;
% Omega
Coef  = {1/dt, k};
Test  = {'v.val', 'v.grad'};
Trial = {'u.val', 'u.grad'};
kk = assem2d(Th,Coef,Test,Trial,Vh,quadOrder); % stiffness matrix
% Robin
Th.elem1d = Th.bdEdgeType{1};
%Th.elem1dIdx = Th.bdEdgeIdxType{1};

Coef = alpha;
Test = 'v.val';
Trail = 'u.val';
kk = kk + assem1d(Th,Coef,Test,Trial,Vh,quadOrder);

%% Time iterations
t = 0;
uh0 = interp2d(u0,Th,Vh);
while t<T

    % Linear form
    Coef = uh0/dt;  Test = 'v.val';
    ff = assem2d(Th,Coef,Test,[],Vh,quadOrder);

    % Neumann boundary condition
    Coef = alpha*ue;  Test = 'v.val';
    ff = ff + assem1d(Th,Coef,Test,[],Vh,quadOrder);
    
    % Dirichlet boundary condition
    on = 2; gD = u0;
    uh = apply2d(on,Th,kk,ff,gD,Vh);

    % Update
    uh0 = uh;
    t = t + dt;

    % Interpolate to a 2D cartesian grid
    figure(1);
    x = 0:0.01:6; y = 0:0.01:1;
    varcontourf(x,y,node,elem,uh(1:Th.N),20);
    color = jet(10);  colormap(color);
    drawnow;
end

%% Show solution
figure, 
subplot(1,2,1), 
showsolution(node,elem,uh);
title('Numerical solution of varFEM');
uff = solFreeFEM('sol_thermalConduction.txt');
subplot(1,2,2),
showsolution(node,elem,uff);
title('Numerical solution of FreeFEM++');