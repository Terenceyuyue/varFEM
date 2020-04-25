function [Ndof,NNdof] = dofnum(Th,Vh)


%% P1-Lagrange
if strcmpi(Vh, 'P1') 
    Ndof = 3; NNdof = Th.N;
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')    
    Ndof = 6; NNdof = Th.N + Th.NE;
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3') 
    Ndof = 10; NNdof = Th.N + 2*Th.NE + Th.NT;
end