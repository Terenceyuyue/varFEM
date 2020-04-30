function u = Applyboundary2D(Th,kk,ff,g_D,Vh)
%Ex1: g_D = pde.g_D;  % scalar
%Ex2: g_D = { g_D1, g_D2, g_D3 }
%Ex3: g_D = { [], [], g_D3 }

%% ------------------ input -----------------------
if nargin==4, Vh = {'P1'}; end  % default: P1
if ~iscell(Vh), Vh = {Vh}; end
nSpace = length(Vh);

%% ------------------ scalar case -----------------------
if nSpace==1
    [bddof,bdval] = getbd(Th,g_D,Vh);
    u = zeros(size(ff)); u(bddof) = bdval;
    ff = ff - kk*u;
    freedof = (~bddof);
    u(freedof) = kk(freedof,freedof)\ff(freedof); % direct solver
    return; % otherwise, vectorized FEM
end

%% ------------------ Vectorized FEM -----------------------
% initialization
bddof = cell(nSpace,1); bdval = cell(nSpace,1);  
for i = 1:nSpace
    [~,NNdofu] = dofnum(Th,Vh{i}); % ui
    bddof{i} = false(NNdofu,1); 
end
% modify i-th cell of bddof and bdval
for i = 1:nSpace
    if ~isempty(g_D{i})
        [bddof{i},bdval{i}] = getbd(Th,g_D{i},Vh{i});
    end
end
% empty cell of bdval is automatically deleted when using concatenation
bddof = vertcat(bddof{:});  bdval = vertcat(bdval{:});

u = zeros(size(ff)); u(bddof) = bdval;
ff = ff - kk*u;
freedof = (~bddof);
u(freedof) = kk(freedof,freedof)\ff(freedof); % direct solver