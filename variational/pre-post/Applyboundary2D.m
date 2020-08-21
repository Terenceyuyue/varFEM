function u = Applyboundary2D(Th,A,b,g_D,Vh)
%% APPLYBOUNDARY2D deals with the Dirichlet boundary conditions of 2D problems
%
%Ex1: g_D = pde.g_D;  % scalar
%Ex2: g_D = { g_D1, g_D2, g_D3 } % (u1,u2,u3)
%Ex3: g_D = { [], [], g_D3 }
%
% Copyright (C) Terence Yu.

%% Input
if nargin==4, Vh = {'P1'}; end  % default: P1
if ~iscell(Vh), Vh = {Vh}; end
nSpace = length(Vh);

%% Scalar case
if nSpace == 1 && ~iscell(g_D)
    g_D = {g_D};
end

%% Vector case
% initialization
bddof = cell(nSpace,1); bdval = cell(nSpace,1); NNdofu = zeros(nSpace,1);
for i = 1:nSpace
    NNdofu(i) = dofnum2d(Th,Vh{i}); % ui
    bddof{i} = false(NNdofu(i),1);
end
% modify i-th cell of bddof and bdval
for i = 1:nSpace
    if ~isempty(g_D{i})
        [bddof{i},bdval{i}] = getbd2D(Th,g_D{i},Vh{i});
    end
end
% empty cell of bdval is automatically deleted when using concatenation
NNdof = length(b);
bddof = vertcat(bddof{:});  bdval = vertcat(bdval{:});
freedof = ~bddof;
u = zeros(NNdof,1); u(bddof) = bdval;
b = b - A*u;

%% Solver
u(freedof) = A(freedof,freedof)\b(freedof);

end

%% getbd2D
function [bddof,bdval] = getbd2D(Th,g_D,Vh)
% g_D: Dirichlet function handle

if nargin==2, Vh = 'P1'; end % default: P1

node = Th.node; bdStruct = Th.bdStruct;
%% P1-Lagrange
if  strcmpi(Vh, 'P1')
    NNdof = Th.N;     
    isBdDof = false(NNdof,1);
    fixedNode = bdStruct.bdNodeIdx;
    isBdDof(fixedNode) = true;
    bddof = (isBdDof);
    bdval = g_D(node(fixedNode,:));
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    NNdof = Th.N + Th.NE;
    isBdDof = false(NNdof,1); 
    
    bdStruct = Th.bdStruct;
    bdNodeIdx = bdStruct.bdNodeIdx; 
    bdEdgeIdxD = bdStruct.bdEdgeIdxD;
    bdEdgeD = bdStruct.bdEdgeD;
    
    fixedNode = [bdNodeIdx; bdEdgeIdxD + Th.N];
    isBdDof(fixedNode) = true;
    bddof = isBdDof;  
    zc = (node(bdEdgeD(:,1),:)+node(bdEdgeD(:,2),:))/2;
    node2 = [node(bdNodeIdx,:); zc]; 
    bdval = g_D(node2);
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3')
    N = Th.N; NE = Th.NE; NT = Th.NT; 
    NNdof = N + 2*NE + NT;
    isBdDof = false(NNdof,1); 
    
    bdStruct = Th.bdStruct;
    bdNodeIdx = bdStruct.bdNodeIdx;
    bdEdgeIdxD = bdStruct.bdEdgeIdxD;
    bdEdgeD = bdStruct.bdEdgeD;
    
    fixedNode = [bdNodeIdx; bdEdgeIdxD+N; bdEdgeIdxD+N+NE];
    isBdDof(fixedNode) = true;
    bddof = (isBdDof);
    
    za = (2*node(bdEdgeD(:,1),:) + node(bdEdgeD(:,2),:))/3;
    zb = (node(bdEdgeD(:,1),:) + 2*node(bdEdgeD(:,2),:))/3;
    node2 = [ node(bdNodeIdx,:); za; zb];
    bdval = g_D(node2);    
end

%% Crouzeix-Raviart linear element
if  strcmpi(Vh, 'CR')
    NNdof = Th.NE;     
    isBdDof = false(NNdof,1);
    fixedNode = bdStruct.bdEdgeIdxD;
    bdEdgeD = bdStruct.bdEdgeD;
    isBdDof(fixedNode) = true;
    bddof = (isBdDof);
    zm = (node(bdEdgeD(:,1),:)+node(bdEdgeD(:,2),:))/2;
    bdval = g_D(zm);
end

end

%%
function [NNdof,Ndof] = dofnum2d(Th,Vh)

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
    
    %% Crouzeix-Raviart linear element
    if strcmpi(Vh, 'CR') 
        Ndof = 3; NNdof = Th.NE;
    end

end
