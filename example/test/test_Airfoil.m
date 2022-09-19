clc; clear; close all
%%% This is an exmple given in FreeFem Documentation: Release 4.6
%   Subsection 2.7 - Irrotational Fan Blade Flow and Thermal effects 
%   (see potential.edp)
%
% Airfoil: Consider a wing profile S in a uniform flow. Infinity will be
% represented by a large circle C where the flow is assumed to be of
% uniform velocity. The domain is outside S.
%

%% Parameters
theta = 8*pi/180; 
lift = theta*0.151952/0.0872665; % lift approximation formula
uinfty1 = cos(theta);   uinfty2 = sin(theta);

%% Mesh
Meshid = 1;
switch Meshid
    case 1  % mesh generated by freefem++        
        [node,elem] = getMeshFreeFEM('meshdata_airfoil.msh');
        %load meshdata_airfoil
    case 2  % mesh generated by pdetool        
        [node,elem] = mesh_naca0012();
end
figure, showmesh(node,elem);
naca12 = @(x) 0.17735*sqrt(x) - 0.075597*x - 0.212836*(x.^2) ...
    + 0.17363*(x.^3) - 0.06254*(x.^4);
x = linspace(0,1,1000)';
hold on
plot(x,naca12(x),'r-','linewidth',1);
plot(x,-naca12(x),'r-','linewidth',1);
hold off
axis([-0.5 1.5 -0.5 0.5]);

% mesh info
bdStr = 'x.^2 + y.^2 > 4.5^2'; % 1-C
Th = FeMesh2d(node,elem,bdStr);

%% Bilinear form
Vh = 'P2';  quadOrder = 7;
Coef  = 1;
Test  = 'v.grad';
Trial = 'u.grad';
kk = assem2d(Th,Coef,Test,Trial,Vh,quadOrder); 

%% Linear form
ff = zeros(size(kk,1),1);

%% Dirichlet boundary conditions
on = [1,2];
gBc1 = @(p) uinfty1*p(:,2) - uinfty2*p(:,1); % on 1-C
gBc2 = @(p) -lift + 0*p(:,1);  % on 2-S
uh = apply2d(on,Th,kk,ff,Vh,gBc1,gBc2);

%% Interpolate to a 2D cartesian grid 
tic;
figure, 
subplot(1,2,1),
x = -0.5:0.005:1.5; y = -0.5:0.005:0.5;
varcontourf(x,y,node,elem,uh(1:Th.N),20); 
title('Numerical solution of varFEM');
subplot(1,2,2),
load sol_airfoil.mat ufh
[node,elem] = getMeshFreeFEM('meshdata_airfoil.msh');
varcontourf(x,y,node,elem,ufh(1:size(node,1)),20); 
title('Numerical solution of FreeFEM++');

fprintf('Time for contour plot is: %.4f s \n', toc);