function u = apply2d(on,Th,kk,ff,g_D,Vh)
%%apply2d deals with Dirichlet boundary conditions in 2D
%
%  Ex1: g_D = pde.g_D;  % scalar
%  Ex2: g_D = { g_D1, g_D2, g_D3 } % (u1,u2,u3)
%  Ex3: g_D = { [], [], g_D3 }
%


%% Input
if nargin==4, Vh = {'P1'}; end  % default: P1
if ~iscell(Vh), Vh = {Vh}; end
nSpace = length(Vh);

%% Transform Scalar case to Vector case
if nSpace == 1 && ~iscell(g_D)
    g_D = {g_D};
end

%% Get bddof, and bdval
% initialization
bddof = cell(nSpace,1);   bdval = cell(nSpace,1); 
NNdofu = zeros(nSpace,1);
for i = 1:nSpace
    NNdofu(i) = dofNum2d(Th,Vh{i}); % ui
end
% modify i-th cell of bddof and bdval
for i = 1:nSpace
    if ~isempty(g_D{i})
        [bddof{i},bdval{i}] = getbd2D(on,Th,g_D{i},Vh{i});
        bddof{i} = bddof{i} + (i>=2)*sum(NNdofu(1:i-1));
    end
end
% empty cell of bdval is automatically deleted when using vertcat
NNdof = length(ff);
bddof = vertcat(bddof{:});  
bdval = vertcat(bdval{:});

%% Apply Dirichlet conditions
u = zeros(NNdof,1); u(bddof) = bdval;
ff = ff - kk*u;

%% Set solver
freedof = true(NNdof,1);
freedof(bddof) = false;
u(freedof) = kk(freedof,freedof)\ff(freedof);

end

%% getbd2D
function [bddof,bdval] = getbd2D(on,Th,g_D,Vh)
% g_D: Dirichlet function handle

if nargin==2, Vh = 'P1'; end % default: P1

node = Th.node; 
%% P1-Lagrange
if  strcmpi(Vh, 'P1')  
    fixedNode = Th.bdNodeIdxType{on};
    bddof = fixedNode;
    bdval = g_D(node(fixedNode,:));
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    bdNodeIdxD = Th.bdNodeIdxType{on}; 
    bdEdgeIdxD = Th.bdEdgeIdxType{on};
    bdEdgeD = Th.bdEdgeType{on};
    
    fixedNode = [bdNodeIdxD; bdEdgeIdxD + Th.N];
    bddof = fixedNode;  
    zc = (node(bdEdgeD(:,1),:)+node(bdEdgeD(:,2),:))/2;
    node2 = [node(bdNodeIdxD,:); zc]; 
    bdval = g_D(node2);
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3')
    N = Th.N; NE = Th.NE;     
    bdNodeIdxD = Th.bdNodeIdxType{on};
    bdEdgeIdxD = Th.bdEdgeIdxType{on};
    bdEdgeD = Th.bdEdgeType{on};
    
    fixedNode = [bdNodeIdxD; bdEdgeIdxD+N; bdEdgeIdxD+N+NE];
    bddof = fixedNode;
    
    za = (2*node(bdEdgeD(:,1),:) + node(bdEdgeD(:,2),:))/3;
    zb = (node(bdEdgeD(:,1),:) + 2*node(bdEdgeD(:,2),:))/3;
    node2 = [node(bdNodeIdxD,:); za; zb];
    bdval = g_D(node2);    
end

% %% Crouzeix-Raviart linear element
% if  strcmpi(Vh, 'CR')
%     NNdof = Th.NE;     
%     isBdDof = false(NNdof,1);
%     fixedNode = bdStruct.bdEdgeIdxD;
%     bdEdgeD = bdStruct.bdEdgeD;
%     isBdDof(fixedNode) = true;
%     bddof = (isBdDof);
%     zm = (node(bdEdgeD(:,1),:)+node(bdEdgeD(:,2),:))/2;
%     bdval = g_D(zm);
% end

end