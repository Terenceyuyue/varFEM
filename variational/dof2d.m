function [elem2dof,Ndof,NNdof] = dof2d(Th,Vh)

node = Th.node; elem = Th.elem; 
N = size(node,1); NT = size(elem,1);

%% P1-Lagrange
if strcmpi(Vh, 'P1') 
    % d.o.f. numbers
    Ndof = 3; NNdof = N;
    % local --> global
    elem2dof = elem;
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    % auxstructure
    auxT = Th.auxT; 
    edge = auxT.edge; NE = size(edge,1);
    elem2edge = auxT.elem2edge;
    % d.o.f. numbers
    Ndof = 6; NNdof = N + NE;
    % local --> global
    elem2dof = [elem, elem2edge + N];
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3')
    % auxstructure
    auxT = Th.auxT; 
    edge = auxT.edge; NE = size(edge,1);
    elem2edge = auxT.elem2edge;
    % d.o.f. numbers
    Ndof = 10; NNdof = N + 2*NE + NT;
    % sgnelem
    bdStruct = Th.bdStruct;
    v1 = [2 3 1]; v2 = [3 1 2];  
    bdIndex = bdStruct.bdIndex; E = false(NE,1); E(bdIndex) = 1;
    sgnelem = sign(elem(:,v2)-elem(:,v1));
    sgnbd = E(elem2edge);    sgnelem(sgnbd) = 1;
    sgnelem(sgnelem==-1) = 0;
    elema = elem2edge + N*sgnelem + (N+NE)*(~sgnelem); % 1/3 point
    elemb = elem2edge + (N+NE)*sgnelem + N*(~sgnelem); % 2/3 point
    % local --> global    
    elem2dof = [elem, elema, elemb, (1:NT)'+N+2*NE];
end