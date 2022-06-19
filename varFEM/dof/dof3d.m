function [elem2dof,Ndof,NNdof] = dof3d(Th,Vh)

%% P1-Lagrange
if strcmpi(Vh, 'P1') 
    Ndof = 4; 
    NNdof = Th.N;
    elem2dof = Th.elem3;
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    Ndof = 10; 
    NNdof = Th.N + Th.NE;
    elem2dof = [Th.elem3, Th.elem2edge + Th.N];
end