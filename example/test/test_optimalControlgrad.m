function test_optimalControlgrad
clc;clear;close all
%%% This is an exmple given in FreeFem Documentation: Release 4.6
%   Subsection 2.15 - Optimal Control (see optimcontrol.edp)
%   For this M-file, the gradient is provided.
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
options = optimoptions('fminunc','Display','iter',...
    'SpecifyObjectiveGradient',true); % gradient is provided
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
function [err,derr] = J(z,ud,Th) 
    Ie = @(p) ((p(:,1)-1).^2 + p(:,2).^2 <=4);
    % J
    uh = PDEcon(z,Th);
    fh = Ie(Th.node).*(uh-ud).^2;
    err = integral2d(Th,fh,Th.Vh,Th.quadOrder);
    % dJ: gradient
    derr = zeros(length(z),1);
    for i = 1:length(z)
        dz = [0,0,0];  dz(i) = 1;
        duh = PDEcongrad(z,dz,uh,Th);
        dfh = Ie(Th.node).*(uh-ud).*duh;
        derr(i) = 2*integral2d(Th,dfh,Th.Vh,Th.quadOrder);
    end 
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

%% Problem for computing the gradient
function duh = PDEcongrad(z,dz,uh,Th)     
    
    % Parameters
    Vh = Th.Vh; quadOrder = Th.quadOrder;    
    Ib = @(p) (p(:,1).^2 + p(:,2).^2 < 1.0001);
    Ic = @(p) ( (p(:,1)+3).^2 + p(:,2).^2 < 1.0001);
    Id = @(p) (p(:,1).^2 + (p(:,2)+3).^2 < 1.0001);
    
    
    % Bilinear form
    Coef = @(p) 1 + z(1)*Ib(p) + z(2)*Ic(p) + z(3)*Id(p);
    Test = 'v.grad';
    Trial = 'du.grad';
    kk = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
    
    % Linear form
    cp = @(p) dz(1)*Ib(p) + dz(2)*Ic(p) + dz(3)*Id(p);
    cc = interp2dMat(cp,'.val',Th,Vh,quadOrder);
    uxc = interp2dMat(uh,'.dx',Th,Vh,quadOrder);
    uyc = interp2dMat(uh,'.dy',Th,Vh,quadOrder);
    Coef = {cc.*uxc, cc.*uyc};  
    Test = {'.dx', '.dy'};
    ff = -assem2d(Th,Coef,Test,[],Vh,quadOrder);
    
    % Dirichlet boundary condition
    gD = @(p) 0*p(:,1);
    on = 1;
    duh = apply2d(on,Th,kk,ff,Vh,gD);
end