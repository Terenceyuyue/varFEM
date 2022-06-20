clc;clear;close all
%%% This is an exmple given in FreeFem Documentation: Release 4.6
%   Subsection 2.10 - The System of Stokes for Fluids
%   The Mini element is considered there, which is now replaced by the 
%   Taylor-Hood element since we do not provide the Mini element at present.

Vh = {'P2','P2','P1'}; % v = [v1,v2,q] --> [v1,v2,v3]
quadOrder = 5;

%% Mesh
[node,elem] = getMeshFreeFEM('meshdata_Stokes.msh');
% load meshdata_Stokes
Th = FeMesh2d(node,elem); 

%% PDE data
pde = Stokesdatavar;

%% Assemble stiffness matrix
vstr = {'v1','v2','q'}; ustr = {'u1','u2','p'};
% [v1,v2,q] = [v1,v2,v3], [u1,u2,p] = [u1,u2,u3]
%   a1(v1,u1) + a2(v2,u2) 
% + b1(v1,p)  + b2(v2,p) 
% + b1(q,u1)  + b2(q,u2) - eps*(q,p)
eps = 1e-10;
Coef = { 1, 1,  -1,-1,  -1,-1,  -eps};  
Test  = {'v1.grad', 'v2.grad', 'v1.dx',  'v2.dy',  'q.val', 'q.val', 'q.val'};
Trial = {'u1.grad', 'u2.grad', 'p.val',  'p.val', 'u1.dx',  'u2.dy',  'p.val'};
[Test,Trial] = getStdvarForm(vstr, Test, ustr, Trial); % [u1,u2,p] --> [u1,u2,u3]
[kk,info] = int2d(Th,Coef,Test,Trial,Vh,quadOrder);

%% Assemble right hand side
trf = eye(2); 
Coef1 = @(pz) pde.f(pz)*trf(:, 1);  Coef2 = @(pz) pde.f(pz)*trf(:, 2);
Coef = {Coef1, Coef2};
Test = {'v1.val', 'v2.val'};
ff = int2d(Th,Coef,Test,[],Vh,quadOrder);

%% Apply Dirichlet boundary conditions
tru = eye(2);
g_D1 = @(pz) pde.g_D(pz)*tru(:, 1);
g_D2 = @(pz) pde.g_D(pz)*tru(:, 2);
g_D = {g_D1, g_D2, []};
on = 1;
U = apply2d(on,Th,kk,ff,Vh,g_D);
NNdofu = info.NNdofu; 
id1 =  NNdofu(1);  id2 = NNdofu(1)+ NNdofu(2); 
uh = [ U(1:id1), U(id1+1:id2) ]; % uh = [u1h, u2h]
ph = U(id2+1:end);

%% Interpolate to a 2D cartesian grid
figure, 
x = 0:0.01:1; y = x;
subplot(1,3,1),
varcontourf(x,y,node,elem,uh(1:Th.N,1),20); title('u_1');
subplot(1,3,2),
varcontourf(x,y,node,elem,uh(1:Th.N,2),20); title('u_2');
subplot(1,3,3),
varcontourf(x,y,node,elem,ph(1:Th.N),20); title('p');