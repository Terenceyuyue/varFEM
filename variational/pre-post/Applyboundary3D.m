function u = Applyboundary3D(Th,A,b,g_D,Vh)
%% APPLYBOUNDARY2D deals with the Dirichlet boundary conditions of 2D problems
%
%Ex1: g_D = pde.g_D;  % scalar
%Ex2: g_D = { g_D1, g_D2, g_D3 } % (u1,u2,u3)
%Ex3: g_D = { [], [], g_D3 }
%
% Copyright (C) Terence Yu.

%% input
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
    [~,NNdofu(i)] = dofnum3D(Th,Vh{i}); % ui
    bddof{i} = false(NNdofu(i),1);
end
% modify i-th cell of bddof and bdval
for i = 1:nSpace
    if ~isempty(g_D{i})
        [bddof{i},bdval{i}] = getbd3D(Th,g_D{i},Vh{i});
    end
end
% empty cell of bdval is automatically deleted when using concatenation
NNdof = length(b);
bddof = vertcat(bddof{:});  bdval = vertcat(bdval{:});
u = zeros(NNdof,1);  u(bddof) = bdval;
b = b - A*u;
freedof = (~bddof);
%% Solver
solver = 'amg';
if NNdof < 2e-3 || nSpace > 1, solver = 'direct'; end
switch solver
    case 'direct'        
        u(freedof) = A(freedof,freedof)\ff(freedof); % direct solver
    case 'amg'
        option.solver = 'CG';
        u(freedof) = amg(A(freedof,freedof),b(freedof),option);  
end
u(bddof) = bdval;

end

%% dofnum3D
function [Ndof,NNdof] = dofnum3D(Th,Vh)

%% P1-Lagrange
if strcmpi(Vh, 'P1')
    Ndof = 4; NNdof = Th.N;
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    Ndof = 10; NNdof = Th.N + Th.NE + Th.NT;
end

end

%% getbd3D
function [bddof,bdval] = getbd3D(Th,g_D,Vh)
% g_D: Dirichlet function handle

if nargin==2, Vh = 'P1'; end % default: P1
node = Th.node;  bdStruct = Th.bdStruct;

%% P1-Lagrange
if  strcmpi(Vh, 'P1')
    NNdof = Th.N;
    isBdDof = false(NNdof,1);
    fixedNode = bdStruct.bdNodeDIdx;
    isBdDof(fixedNode) = true;
    bddof = (isBdDof);
    bdval = g_D(node(fixedNode,:));
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    N = Th.N;  NE = Th.NE;  NNdof = N + NE;    
    isBdDof = false(NNdof,1);
    
    bdNodeDIdx = bdStruct.bdNodeDIdx;
    bdEdgeDIdx = bdStruct.bdEdgeDIdx;
    bdEdgeD = bdStruct.bdEdgeD;
    fixedNode = [bdNodeDIdx; bdEdgeDIdx+N];
    isBdDof(fixedNode) = true;
    bddof = (isBdDof);
    
    node2 = [node(bdNodeDIdx,:);
        (node(bdEdgeD(:,1),:) + node(bdEdgeD(:,2),:))/2];
    bdval = g_D(node2);
    
end
end

