function [elem2dof,ndof,NNdof] = dof1d(Th,feSpace)

%% mesh information
node = Th.node;  N = size(node,1);
elem1D = Th.elem1D; nel = size(elem1D,1);
if isfield(Th, 'bdIndex1D'), bdIndex1D = Th.bdIndex1D; end

%% P1-Lagrange
if  strcmpi(feSpace, 'P1')
    % d.o.f. numbers
    ndof = 2; NNdof = N;
    % local --> global
    elem2dof = elem1D;
end

%% P2-Lagrange
if  strcmpi(feSpace, 'P2')
    if isfield(Th, 'elem') % 2D
        % auxstructure
        auxT = Th.auxT; edge = auxT.edge;
        NE = size(edge,1);
        % d.o.f. numbers
        ndof = 3; NNdof = N + NE;
        % local --> global
        elem2dof = [elem1D, bdIndex1D + N];
    else
        ndof = 3; NNdof = N + nel;
        elem2dof = [elem1D, (1:nel)'+N];
    end
end

%% P3-Lagrange
if  strcmpi(feSpace, 'P3')
    if isfield(Th, 'elem') % 2D
        % auxstructure
        auxT = Th.auxT; edge = auxT.edge;
        NT = size(Th.elem,1);  NE = size(edge,1);
        % d.o.f. numbers
        ndof = 4; NNdof = N + 2*NE + NT;
        % local --> global
        elem2dof = [elem1D, bdIndex1D+N, bdIndex1D+N+NE];
    else
        ndof = 4; NNdof = N + 2*nel;
        elem2dof = [elem1D, (1:nel)'+N, (1:nel)'+N+nel];
    end
end