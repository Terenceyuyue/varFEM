function [elem2dof,Ndof,NNdof] = dof2d(Th,Vh)
%% DOF2D returns elem2dof for assembling 2D part

%% Mesh information
bdType = 2;
if isfield(Th,'bdType'), bdType = Th.bdType; end
node = Th.node;  elem = Th.elem;

%% P1-Lagrange
if strcmpi(Vh, 'P1')
    if size(node,2) == 2 % 2-D
        % d.o.f. numbers
        N = Th.N;
        Ndof = 3; NNdof = N;
        % local --> global
        elem2dof = elem;
    end
    if size(node,2) == 3 % 3-D boundary conditions
        % elem2dof in 3-D
        elem2dof3 = dof3d(Th,Vh);
        % d.o.f. numbers
        N = Th.N;
        Ndof = 3; NNdof = N;
        % elem2dof in 2-D inherited from 3-D
        bdFlag = Th.bdFlag;
        elem2dof = [elem2dof3(bdFlag(:,1)==bdType, [2,3,4]);
            elem2dof3(bdFlag(:,2)==bdType, [1,4,3]);
            elem2dof3(bdFlag(:,3)==bdType, [1,2,4]);
            elem2dof3(bdFlag(:,4)==bdType, [1,3,2])];
    end
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    if size(node,2) == 2 % 2-D
        % auxstructure
        elem2edge = Th.elem2edge;
        % d.o.f. numbers
        N = Th.N; NE = Th.NE;
        Ndof = 6; NNdof = N + NE;
        % local --> global
        elem2dof = [elem, elem2edge + N];
    end
    if size(node,2) == 3 % 3-D boundary conditions
        % elem2dof in 3-D
        elem2dof3 = dof3d(Th,Vh);
        % d.o.f. numbers
        N = Th.N; NE = Th.NE;
        Ndof = 6; NNdof = N + NE;
        % elem2dof in 2-D inherited from 3-D
        bdFlag = Th.bdFlag;
        elem2dof = [elem2dof3(bdFlag(:,1)==bdType, [2,3,4,10,9,8]);
            elem2dof3(bdFlag(:,2)==bdType, [1,4,3,10,6,7]);
            elem2dof3(bdFlag(:,3)==bdType, [1,2,4,9,7,5]);
            elem2dof3(bdFlag(:,4)==bdType, [1,3,2,8,5,6])];
    end
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3')
    if size(node,2) == 2 % 2-D
        % auxstructure
        N = Th.N; NT = Th.NT; NE = Th.NE;
        elem2edge = Th.elem2edge;
        % d.o.f. numbers
        Ndof = 10; NNdof = N + 2*NE + NT;
        % sgnelem
        bdStruct = Th.bdStruct;
        v1 = [2 3 1]; v2 = [3 1 2];
        bdEdgeIdx = bdStruct.bdEdgeIdx; E = false(NE,1); E(bdEdgeIdx) = 1;
        sgnelem = sign(elem(:,v2)-elem(:,v1));
        sgnbd = E(elem2edge);    sgnelem(sgnbd) = 1;
        sgnelem(sgnelem==-1) = 0;
        elema = elem2edge + N*sgnelem + (N+NE)*(~sgnelem); % 1/3 point
        elemb = elem2edge + (N+NE)*sgnelem + N*(~sgnelem); % 2/3 point
        % local --> global
        elem2dof = [elem, elema, elemb, (1:NT)'+N+2*NE];
    end
    % P3-Lagrange is not provided in 3-D
end