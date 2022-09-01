function [NNdof,Ndof] = dofNum3d(Th,Vh)

%% P1-Lagrange
if strcmpi(Vh, 'P1')
    Ndof = 4; 
    NNdof = Th.N;
    return;
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    Ndof = 10; 
    NNdof = Th.N + Th.NE;
    return;
end