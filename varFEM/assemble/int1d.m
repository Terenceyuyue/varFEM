function [output,varargout] = int1d(Th,Coef,Test,Trial,Vh,quadOrder)
%%int1d assembles bilinear or linear forms of scalar or vector type 
% e.g.
% Vh --> (v1,v2,v3) or (u1,u2,u3) --> { 'P2','P2','P1' }
% vi and ui are in the same FE spaces for i = 1,2,3
%

%% Default input for FE space and quadrature order
if nargin == 3,  Trial = []; Vh = {'P1'}; quadOrder = 5; end  % linear part: Trial = []
if nargin == 4,  Vh = {'P1'}; quadOrder = 5; end % default: P1
if nargin == 5, quadOrder = 5; end

if ~iscell(Vh), Vh = {Vh}; end % Vh is stored in cell array
if isnumeric(Coef) && length(Coef)==1 % a constant
    Coef = @(p) Coef+0*p(:,1);
end
if ~iscell(Coef), Coef = {Coef}; end 
if ~isempty(Test) && ~iscell(Test), Test = {Test}; end
if ~isempty(Trial) && ~iscell(Trial), Trial = {Trial}; end

%% Scalar case
nSpace = length(Vh);
if nSpace == 1
    output = assem1d(Th,Coef,Test,Trial,Vh,quadOrder);
    return;  % otherwise, vector case
end

%% extended [Coef,Test,Trial]
% e.g. Eij(u)Eij(v) = v1x*u1x + v2x*u2y + 0.5(v1y+v2x)*(u1y+u2x)
% is expanded as
%  Coef  = {1, 1, 0.5, 0.5, 0.5, 0.5};
%  Test  = {'v1.dx', 'v2.dy', 'v1.dy', 'v1.dy', 'v2.dx', 'v2.dx'};
%  Trial = {'u1.dx', 'u2.dy', 'u1.dy', 'u2.dx', 'u1.dy', 'u2.dx'};
% nSpace = length(Vh);   % nSpace = 3 for (v1,v2,v3)
if ~isempty(Trial) && nSpace>1 && sum(mycontains(Trial,'+'))
    [Coef,Test,Trial] = getExtendedvarForm(Coef,Test,Trial);
end

%% Sparse assembly index w.r.t. the vectorized FE space
% elementwise d.o.f.s
elem2dofv = cell(1,nSpace); NNdofv = zeros(1,nSpace);
for i = 1:nSpace
    [elem2dofv{i},~,NNdofv(i)] = dof1d(Th,Vh{i}); % vi
end
elem2dofu = elem2dofv; NNdofu = NNdofv;
NNdofvv = sum(NNdofv);  NNdofuu = NNdofvv;

info.NNdofu = NNdofu; 
varargout{1} = info;

% assembly index
ii = cell(nSpace^2,1);  jj = ii;  s = 1;
for i = 1:nSpace
    for j = 1:nSpace
        [iiv,jju] = getSparse(elem2dofv{i}, elem2dofu{j});
        ii{s} = iiv + (i>=2)*sum(NNdofv(1:i-1));
        jj{s} = jju + (j>=2)*sum(NNdofu(1:j-1));
        s = s+1;
    end
end

%% Bilinear form
if ~isempty(Trial)
    idvu = true(nSpace^2,1);  ss = cell(nSpace^2,1); k = 1;
    % (vi,uj), i,j = 1,...,nSpace
    for i = 1:nSpace
        for j = 1:nSpace
            Vhij = { Vh{i}, Vh{j} }; % FE space pair w.r.t. (vi,uj)
            id = mycontains(Test,num2str(i)) & mycontains(Trial,num2str(j));
            % id gives all the pairs having (vi,uj)
            if sum(id)==0  %  empty
                idvu(k) = false; k = k+1;   continue;
            end
            [~, ss{k}] = assem1d(Th,Coef(id),Test(id),Trial(id),Vhij,quadOrder);
            k = k+1;
        end
    end
    ii = ii(idvu); jj = jj(idvu); ss = ss(idvu);
    ii = vertcat(ii{:});   jj = vertcat(jj{:});  ss = vertcat(ss{:});
    output = sparse(ii,jj,ss,NNdofvv,NNdofuu);
    return;  % The remaining code will be neglected.
end

%% Linear form
ff = zeros(NNdofvv,1);

% Test --> v.val = [v1.val, v2.val, v3.val]
if length(Test)==1 && strcmpi(Test{1}, 'v.val')
    if length(Coef) == 1  % pde.f
        trf = eye(nSpace); f = Coef{1};
        Coef = cell(nSpace,1);
        for i = 1:nSpace
            Coef{i} = @(pz) f(pz)*trf(:,i); 
        end
    end  % else: Coef = {f1Mat,f2,f3}
    for i = 1:nSpace
        Coefv = Coef{i};
        Testv = sprintf('v%d.val',i);
        Fi = assem1d(Th,Coefv,Testv,[],Vh{i},quadOrder);
        id = (1:NNdofv(i)) + (i>=2)*sum(NNdofv(1:i-1)); 
        ff(id) = ff(id) + Fi;
    end
    output = ff;  return;
end

% Test -->  'v1.val' --> {'v1.val'}
if ~iscell(Test), Coef = {Coef};  Test = {Test}; end

% Test --> {'v1.val', 'v3.val'}
for s = 1:length(Test)
    Coefv = Coef{s}; 
    Testv = Test{s};
    for i = 1:nSpace
        str = sprintf('v%d',i);  % vi
        if mycontains(Testv, str)
            Fi = assem1d(Th,Coefv,Testv,[],Vh{i},quadOrder);
            id = (1:NNdofv(i)) + (i>=2)*sum(NNdofv(1:i-1)); 
            ff(id) = ff(id) + Fi;
        end
    end
end
output = ff;