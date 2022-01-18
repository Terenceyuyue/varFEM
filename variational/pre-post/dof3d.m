function [elem2dof,Ndof,NNdof] = dof3d(Th,Vh)
%% DOF3D returns elem2dof for assembling 3D part
%
% Copyright (C) Terence Yu.

%% P1-Lagrange
if strcmpi(Vh, 'P1') 
    % d.o.f. numbers
    Ndof = 4; NNdof = Th.N;
    % local --> global
    elem2dof = Th.elem;
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    % d.o.f. numbers
    Ndof = 10; NNdof = Th.N + Th.NE;
    % local --> global
    elem2dof = [Th.elem, Th.elem2edge + Th.N];
end