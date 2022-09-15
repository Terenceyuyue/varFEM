function [kkDof,freedof] = apply2dMat(kk,on,Th,Vh,varargin)

if ~iscell(Vh), Vh = {Vh}; end

%% Input
nBc = length(on);
if nargin==1
    gDLogic = varargin{1};
    gBc = repmat({gDLogic},1,nBc);
else
    gBc = varargin;
end
nSpace = length(Vh);

%% Get bddof, and bdval
% initialization
bddof = cell(nBc,nSpace);   
NNdofu = zeros(nSpace,1);
for i = 1:nSpace
    NNdofu(i) = dofNum2d(Th,Vh{i}); % ui
end
% modify i-th cell of bddof and bdval
for s = 1:nBc  % on = s
    gDLogic = gBc{s};
    if ~iscell(gDLogic), gDLogic = {gDLogic}; end
    for i = 1:nSpace  % g_D = {g_D1, [], g_D3}
        if ~isempty(gDLogic{i})
            bddof{s,i} = getbd2D(on(s),Th,Vh{i});
            bddof{s,i} = bddof{s,i} + (i>=2)*sum(NNdofu(1:i-1));
        end
    end
end
% empty cell of bdval is automatically deleted when using vertcat
NNdof = size(kk,1);
bddof = vertcat(bddof{:});  

%% The new matrix
freedof = true(NNdof,1);  freedof(bddof) = false;
kkDof = kk(freedof,freedof);

end

%% getbd2D
function bddof = getbd2D(on,Th,Vh)
% g_D: Dirichlet function handle

if nargin==2, Vh = 'P1'; end % default: P1

%% P1-Lagrange
if  strcmpi(Vh, 'P1')  
    fixedNode = Th.bdNodeIdxType{on};
    bddof = fixedNode;
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    bdNodeIdxD = Th.bdNodeIdxType{on}; 
    bdEdgeIdxD = Th.bdEdgeIdxType{on};
    
    fixedNode = [bdNodeIdxD; bdEdgeIdxD + Th.N];
    bddof = fixedNode;  
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3')
    N = Th.N; NE = Th.NE;     
    bdNodeIdxD = Th.bdNodeIdxType{on};
    bdEdgeIdxD = Th.bdEdgeIdxType{on};
    
    fixedNode = [bdNodeIdxD; bdEdgeIdxD+N; bdEdgeIdxD+N+NE];
    bddof = fixedNode;
end

end