function test_optimalControl
clc;clear;close all
%%% This is an exmple given in FreeFem Documentation: Release 4.6
%   Subsection 2.15 - Optimal Control (see optimcontrol.edp)
%   For this M-file, the gradient is not provided.
%

%% Mesh
[node,elem] = getMeshFreeFEM('meshdata_optimcontrol.msh');
%load meshdata_optimcontrol node elem
%showmesh(node,elem);
% mesh info
Th = FeMesh2d(node,elem);
% Vh, quadOrder
Th.Vh = 'P1';  Th.quadOrder= 5;

%% The constructed solution
zd = [2, 3, 4];
ud = PDEcon(zd,Th);

%% Find the mimimizer
options.LargeScale = 'off';
options.HessUpdate = 'bfgs';
options.Display = 'iter';
z0 = [1,1,1];
zmin = fminunc(@(z) J(z,ud,Th), z0, options)

%% Plot
uh = PDEcon(zmin,Th);
uff = solFreeFEM('sol_uoptimcontrol.txt');   % zd = [2, 3, 4]
figure,
subplot(2,2,1),
showsolution(node,elem,uh);
subplot(2,2,3)
showsolution(node,elem,uff);
subplot(2,2,2),
x = -5:0.05:5; y = x;
varcontourf(x,y,node,elem,uh,20);
title('Numerical solution of varFEM');
subplot(2,2,4),
varcontourf(x,y,node,elem,uff(1:Th.N),20);
title('Numerical solution of FreeFEM++');

end

%% Cost function
function err = J(z,ud,Th) 
    Ie = @(p) ((p(:,1)-1).^2 + p(:,2).^2 <=4);
    uh = PDEcon(z,Th);

    fh = Ie(Th.node).*(uh-ud).^2;
    err = integral2d(Th,fh,Th.Vh,Th.quadOrder);
end

%% Problem for the PDE constraint
function uh = PDEcon(z,Th)     
    
    % Parameters
    Vh = Th.Vh; quadOrder = Th.quadOrder;    
    Ib = @(p) (p(:,1).^2 + p(:,2).^2 < 1.0001);
    Ic = @(p) ( (p(:,1)+3).^2 + p(:,2).^2 < 1.0001);
    Id = @(p) (p(:,1).^2 + (p(:,2)+3).^2 < 1.0001);
    
    
    % Bilinear form
    Coef = @(p) 1 + z(1)*Ib(p) + z(2)*Ic(p) + z(3)*Id(p);
    Test = 'v.grad';
    Trial = 'u.grad';
    kk = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
    
    % Linear form
    ff = zeros(size(kk,1),1);
    
    % Dirichlet boundary condition
    gD = @(p) p(:,1).^3 - p(:,2).^3;
    on = 1;
    uh = apply2d(on,Th,kk,ff,Vh,gD);
end