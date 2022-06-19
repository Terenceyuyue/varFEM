function [elem2dof,Ndof,NNdof] = dof2d(Th,Vh)

%% Mesh info
node = Th.node;
if size(node,2)==2
    elem = Th.elem;
end
if size(node,2)==3 && isfield(Th,'on')
    on = Th.on;
    elem = Th.bdFaceType{on};
    bdFace2edge = Th.bdFace2edgeType{on};
end
if isfield(Th,'elem'), elem = Th.elem; end % given by userTh.on = 1;  
if isfield(Th,'bdFace2edge'), bdFace2edge = Th.bdFace2edge; end % given by userTh.on = 1;  

%% P1-Lagrange
if strcmpi(Vh, 'P1')
    if size(node,2) == 2 % 2-D
        Ndof = 3; % number of local dofs
        NNdof = size(node,1); % number of global dofs
        elem2dof = elem; 
    end
    if size(node,2) == 3 % 3-D boundary faces
        Ndof = 3;
        NNdof = Th.N;
        elem2dof = elem; % elem2dof in 2-D inherited from 2-D
    end
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    if size(node,2) == 2 % 2-D
        Ndof = 6;
        NNdof = Th.N + Th.NE;
        elem2dof = [elem, Th.elem2edge+Th.N];
    end
    if size(node,2)==3
        Ndof = 6; 
        NNdof = Th.N + Th.NE;
        elem2dof = [elem, bdFace2edge + Th.N];
    end
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3')
    if size(node,2) == 2 % 2-D
        % dof numbers
        N = Th.N; NT = Th.NT; NE = Th.NE;
        Ndof = 10; NNdof = N + 2*NE + NT;
        % sgnelem
        elem2edge = Th.elem2edge;  bdEdgeIdx = Th.bdEdgeIdx; 
        v1 = [2 3 1]; v2 = [3 1 2];        
        E = false(NE,1); E(bdEdgeIdx) = 1;
        sgnelem = sign(elem(:,v2)-elem(:,v1));
        sgnbd = E(elem2edge);    sgnelem(sgnbd) = 1;
        sgnelem(sgnelem==-1) = 0;
        elema = elem2edge + N*sgnelem + (N+NE)*(~sgnelem); % 1/3 point
        elemb = elem2edge + (N+NE)*sgnelem + N*(~sgnelem); % 2/3 point
        % elem2dof
        elem2dof = [elem, elema, elemb, (1:NT)'+N+2*NE];
    end
    % P3-Lagrange is not provided in 3-D
end