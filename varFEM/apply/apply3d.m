function u = apply3d(on,Th,kk,ff,Vh,varargin)
%%apply2d deals with Dirichlet boundary conditions in 2D
%
%  varargin = gBc1, gBc2 if on = [1,2]
%
%  gBc1 or gBc2 has the following forms:
%    Ex1: g_D = pde.g_D;  % scalar
%    Ex2: g_D = { g_D1, g_D2, g_D3 } % (u1,u2,u3)
%    Ex3: g_D = { [], [], g_D3 }
%

if ~iscell(Vh), Vh = {Vh}; end

%% Input
nBc = length(on);
if nargin==1
    g_D = varargin{1};
    gBc = repmat({g_D},1,nBc);
else
    gBc = varargin;
end
nSpace = length(Vh);

%% Get bddof, and bdval
% initialization
bddof = cell(nBc,nSpace);   bdval = cell(nBc,nSpace); 
NNdofu = zeros(nSpace,1);
for i = 1:nSpace
    NNdofu(i) = dofNum3d(Th,Vh{i}); % ui
end
% modify i-th cell of bddof and bdval
for s = 1:nBc  % on = s
    g_D = gBc{s};
    if ~iscell(g_D), g_D = {g_D}; end
    for i = 1:nSpace  % g_D = {g_D1, [], g_D3}
        if ~isempty(g_D{i})
            [bddof{s,i},bdval{s,i}] = getbd3D(on(s),Th,g_D{i},Vh{i});
            bddof{s,i} = bddof{s,i} + (i>=2)*sum(NNdofu(1:i-1));
        end
    end
end
% empty cell of bdval is automatically deleted when using vertcat
NNdof = length(ff);
bddof = vertcat(bddof{:});  
bdval = vertcat(bdval{:});

%% Apply Dirichlet conditions
u = zeros(NNdof,1); u(bddof) = bdval;
ff = ff - kk*u;

%% Solver
freedof = true(NNdof,1);  freedof(bddof) = false;
if ~isfield(Th,'solver')
    solver = 'amg';
else
    solver = Th.solver;
end
if NNdof < 2e3, solver = 'direct'; end
switch solver
    case 'direct'        
        u(freedof) = kk(freedof,freedof)\ff(freedof); % direct solver
    case 'amg'
        option.solver = 'CG';
        u(freedof) = amg(kk(freedof,freedof),ff(freedof),option);  
    case 'cg'
        u(freedof) = cgs(kk(freedof,freedof),ff(freedof),1e-12,NNdof);
end
u(bddof) = bdval;

end


%% getbd3D
function [bddof,bdval] = getbd3D(on,Th,g_D,Vh)
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
    N = Th.N;  
    bdNodeIdxD = Th.bdNodeIdxType{on};
    bdEdgeIdxD = Th.bdEdgeIdxType{on};
    bdEdgeD = Th.bdEdgeType{on};

    fixedNode = [bdNodeIdxD; bdEdgeIdxD+N];
    bddof = fixedNode;
    
    node2 = [node(bdNodeIdxD,:);
        (node(bdEdgeD(:,1),:) + node(bdEdgeD(:,2),:))/2];
    bdval = g_D(node2);    
end

end

