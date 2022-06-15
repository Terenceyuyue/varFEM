clc; clear; close all
%%% This is an exmple given in FreeFEM Documentation: Release 4.6
%   Subsection 2.3 - Membrane (membrane.edp)

%% Mesh
% ellipse with a = 2, b = 1
a = 2; b = 1; 
g = ellipseg(a,b);
[p,e,t] = initmesh(g,'hmax',0.2);
[p,e,t] = refinemesh(g,p,e,t);
node = p'; elem = t(1:3,:)';
figure, 
subplot(1,2,1),
showmesh(node,elem);
% bdStr
bdNeumann = 'y<0 & x>-sin(pi/3)'; % string for Neumann
% mesh info
Th = FeMesh2d(node,elem,bdNeumann);

%% PDE data
f = 1; % right-hand side
g_D = @(p) p(:,1); % Dirichlet: gD = x  
g_N = 0; % Neumann: no need to assemble it

%% Assemble bilinear form
Vh = 'P1'; quadOrder = 5;
Coef = 1;  Test = 'v.grad';  Trial = 'u.grad';
kk = assem2d(Th,Coef,Test,Trial,Vh,quadOrder); 

%% Assemble the right-hand side
Coef = f;  Test = 'v.val';
ff = assem2d(Th,Coef,Test,[],Vh,quadOrder);

%% Apply Dirichlet boundary value conditions
on = 2;
uh = apply2d(on,Th,kk,ff,Vh,g_D);

subplot(1,2,2),
showsolution(node,elem,uh);

%% Interpolate to a 2D cartesian grid
figure, 
varcontourf(node,elem,uh,20);