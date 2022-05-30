function [elem2dof,Ndof,NNdof] = dof1d(Th,Vh)
%% DOF1D returns elem2dof for assembling 1D part
%

% Copyright (C) Terence Yu.


node = Th.node;  elem1d = Th.elem1d;
N = size(node,1); nel = size(elem1d,1);

%% P1-Lagrange
if  strcmpi(Vh, 'P1')
    if size(node,2) == 1 % 1-D
        Ndof = 2;
        NNdof = N; % number of global dofs
        elem2dof = elem1d;        
    end
    if size(Th.node,2) == 2 % 2-D boundary condition
        Ndof = 2; 
        NNdof = N;        
        elem2dof = elem1d;  % elem2dof in 1-D inherited from 2-D
    end    
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    if size(node,2) == 1  % 1-D
        Ndof = 3; 
        NNdof = N + nel;
        elem2dof = [elem1d, (1:nel)'+N];
    end
    if size(node,2) == 2  % 2-D boundary condition
        NE = Th.NE;
        Ndof = 3; 
        NNdof = N + NE;
        elem2dof = [elem1d, Th.bdEdgeIdx1d + N];
    end    
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3')
    if size(Th.node,2) == 1  % 1-D
        Ndof = 4; 
        NNdof = N + 2*nel;
        elem2dof = [elem1d, (1:nel)'+N, (1:nel)'+N+nel];
    end
    if size(Th.node,2) == 2  % 2-D boundary condition
        % dof numbers
        N = Th.N; NE = Th.NE; NT = Th.NT;
        Ndof = 4; 
        NNdof = N + 2*NE + NT;
        % elem2dof
        elem1d = Th.elem1d;
        bdEdgeIdx1d = Th.bdEdgeIdx1d;
        elem2dof = [elem1d, bdEdgeIdx1d+N, bdEdgeIdx1d+N+NE];
    end
end