function [elem2dof,Ndof,NNdof] = dof1d(Th,Vh)
%% DOF1D returns elem2dof for assembling 1D part

%% Mesh information
bdType = 2;
if isfield(Th,'bdType'), bdType = Th.bdType; end
elem = Th.elem;
if size(Th.node,2) == 1, nel = size(elem,1); end % 1-D

%% P1-Lagrange
if  strcmpi(Vh, 'P1')
    if size(Th.node,2) == 1 % 1-D
        N = Th.N;
        % d.o.f. numbers
        Ndof = 2; NNdof = N;
        % local --> global
        elem2dof = elem;        
    end
    if size(Th.node,2) == 2 % 2-D boundary conditions
        % d.o.f. numbers
        N = Th.N; 
        Ndof = 2; NNdof = N;
        % elem2dof in 1-D inherited from 2-D
        bdStruct = Th.bdStruct;
        elem1D = bdStruct.bdEdgeN;
        if bdType==1
            elem1D = bdStruct.bdEdgeD;
        end
        elem2dof = elem1D;
    end    
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    if size(Th.node,2) == 1  % 1-D
        N = size(Th.node,1);
        Ndof = 3; NNdof = N + nel;
        elem2dof = [elem, (1:nel)'+N];
    end
    if size(Th.node,2) == 2  % 2-D boundary conditions
        % d.o.f. numbers
        N = Th.N; NE = Th.NE;
        Ndof = 3; NNdof = N + NE;
        % local --> global
        bdStruct = Th.bdStruct;
        elem1D = bdStruct.bdEdgeN;
        bdEdgeIdx1D = bdStruct.bdEdgeIdxN;
        if bdType==1
            elem1D = bdStruct.bdEdgeD;
            bdEdgeIdx1D = bdStruct.bdEdgeIdxD;
        end
        elem2dof = [elem1D, bdEdgeIdx1D + N];
    end
    
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3')
    if size(Th.node,2) == 1  % 1-D
        N = size(Th.node,1);
        Ndof = 4; NNdof = N + 2*nel;
        elem2dof = [elem, (1:nel)'+N, (1:nel)'+N+nel];
    end
    if size(Th.node,2) == 2  % 2-D boundary conditions
        % d.o.f. numbers
        N = Th.N; NE = Th.NE; NT = Th.NT;
        Ndof = 4; NNdof = N + 2*NE + NT;
        % local --> global
        bdStruct = Th.bdStruct;
        elem1D = bdStruct.bdEdgeN;
        bdEdgeIdx1D = bdStruct.bdEdgeIdxN;
        if bdType==1
            elem1D = bdStruct.bdEdgeD;
            bdEdgeIdx1D = bdStruct.bdEdgeIdxD;
        end
        elem2dof = [elem1D, bdEdgeIdx1D+N, bdEdgeIdx1D+N+NE];
    end
end