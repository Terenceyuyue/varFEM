function [elem2dof,Ndof,NNdof] = dof2d(Th,Vh)
%% DOF2D returns elem2dof for assembling 2D part
%
% Copyright (C) Terence Yu.

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
        % d.o.f. numbers
        N = Th.N;
        Ndof = 3; NNdof = N;
        % elem2dof in 2-D inherited from 3-D
        bdStruct = Th.bdStruct;
        if bdType == 1, elem2dof = bdStruct.bdFaceD; end
        if bdType == 2, elem2dof = bdStruct.bdFaceN; end
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
        elem2dof = [elem, elem2edge+N];
    end
    if size(node,2) == 3 % 3-D boundary conditions
        % d.o.f. numbers
        N = Th.N; NE = Th.NE;
        Ndof = 6; NNdof = N + NE;
        % elem2dof in 2-D inherited from 3-D
        bdStruct = Th.bdStruct;
        if bdType == 1
            bdFaceD = bdStruct.bdFaceD;
            bdFace2edgeD = bdStruct.bdFace2edgeD;
            elem2dof = [bdFaceD, bdFace2edgeD+N];
        end
        if bdType == 2
            bdFaceN = bdStruct.bdFaceN;
            bdFace2edgeN = bdStruct.bdFace2edgeN;
            elem2dof = [bdFaceN, bdFace2edgeN+N];
        end
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

%% Crouzeix-Raviart linear element
if strcmpi(Vh, 'CR')
    if size(node,2) == 2 % 2-D
        % d.o.f. numbers
        Ndof = 3; NNdof = Th.NE;
        % local --> global
        elem2dof = Th.elem2edge;
    end
end

%% Morley
if  strcmpi(Vh, 'Morley')
    if size(node,2) == 2 % 2-D
        % auxstructure
        N = Th.N; NE = Th.NE;
        elem2edge = Th.elem2edge;
        % d.o.f. numbers
        Ndof = 6; NNdof = N + NE;
        % local --> global
        elem2dof = [elem, elem2edge+N];
    end
end

%% Zienkiewicz
if  strcmpi(Vh, 'Zienkiewicz')
    if size(node,2) == 2 % 2-D
       N = size(node,1); 
       Ndof = 9; NNdof = 3*N;
       elem2dof = [elem, elem+N, elem+2*N];
    end
end