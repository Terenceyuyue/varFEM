clc;clear;close all

Vh = 'P1'; quadOrder = 5;
Vhvec = repmat( {Vh}, 1, 2 ); % v = [v1,v2]

%% Parameters
E = 21e5; nu = 0.28; 
lambda =  E*nu/((1 + nu)*(1 - 2*nu));
mu = E/(2*(1 + nu));
f = -1;

%% Mesh
[node,elem] = getMeshFreeFEM('meshdata_lame.msh');
% mesh info
bdStr = 'x==0'; % 1-Dirichlet
Th = FeMesh2d(node,elem,bdStr);
% 1-D mesh
Th.on = 2;

%% Assemble stiffness matrix 
% (Eij(u):Eij(v))
Coef = { 1, 1, 0.5 }; 
Test  = {'v1.dx', 'v2.dy', 'v1.dy + v2.dx'};
Trial = {'u1.dx', 'u2.dy', 'u1.dy + u2.dx'};
A = int2d(Th,Coef,Test,Trial,Vhvec,quadOrder);
A = 2*mu*A;

% (div u,div v) 
Coef = { 1 }; 
Test  = { 'v1.dx + v2.dy' };
Trial = { 'u1.dx + u2.dy' };
B = int2d(Th,Coef,Test,Trial,Vhvec,quadOrder);
B = lambda*B;

% stiffness matrix
kk = A + B;

%% Assemble the right hand side 
Coef = f;  Test = 'v2.val';
ff = int2d(Th,Coef,Test,[],Vhvec,quadOrder);

%% Assemble Neumann boundary conditions
% g_N = 0

%% Apply Dirichlet boundary conditions
on = 1;
g_D1 = @(p) 0*p(:,1);
g_D2 = @(p) 0*p(:,1);
g_D = {g_D1, g_D2};
uh = apply2d(on,Th,kk,ff,Vhvec,g_D);

%% Plot
uh = reshape(uh,[],2);
%figure, showsolution(node,elem,uh(:,1));
figure, 
subplot(2,1,1),
options.facecolor = 'w';
showmesh(node,elem,options); hold on
x = linspace(0,20,31); 
y = linspace(-1,1,11);
varquiver(x,y,node,elem,uh,'r','linewidth',2);

%% Movemesh
cf = 100;
node = node + cf*uh;
subplot(2,1,2),
showmesh(node,elem);