function [output,varargout] = int2d(Th,Coef,Test,Trial,Vh,quadOrder)
%% INT2D  assembles bilinear or linear forms of scalar or vector type 
% e.g.
% Vh --> (v1,v2,v3) or (u1,u2,u3) --> { 'P2','P2','P1' }
% vi and ui are in the same FE space
%
% Copyright (C) Terence Yu.

%% Preparation for the input
% default para
if nargin == 3,  Trial = []; Vh = {'P1'}; quadOrder = 3; end
if nargin == 4,  Vh = {'P1'}; quadOrder = 3; end % default: P1
if nargin == 5, quadOrder = 3; end

if ~iscell(Vh), Vh = {Vh}; end % feSpace = 'P1'

%% Scalar case
nSpace = length(Vh);
if nSpace == 1
    output = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
    return;  % otherwise, vecotorial case
end

%% extended [Coef,Trial,Test]
% e.g. Eij(u)Eij(v) = v1x*u1x + v2x*u2y + 0.5(v1y+v2x)*(u1y+u2x)
% is expanded as
%  Coef  = {1, 1, 0.5, 0.5, 0.5, 0.5};
%  Test  = {'v1.dx', 'v2.dy', 'v1.dy', 'v1.dy', 'v2.dx', 'v2.dx'};
%  Trial = {'u1.dx', 'u2.dy', 'u1.dy', 'u2.dx', 'u1.dy', 'u2.dx'};
if ~isempty(Trial) && nSpace>1
    [Coef,Test,Trial] = getExtendedvarForm(Coef,Test,Trial);
end

%% Sparse assembling index of Pk-Lagrange element
% elementwise d.o.f.s
elem2dofv = cell(1,nSpace); NNdofv = zeros(1,nSpace);
for i = 1:nSpace
    [elem2dofv{i},~,NNdofv(i)] = dof2d(Th,Vh{i}); % vi
end
elem2dofu = elem2dofv; NNdofu = NNdofv;
NNdofvv = sum(NNdofv);  NNdofuu = NNdofvv;

info.NNdofu = NNdofu; 
varargout{1} = info;

% assembling index
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
            % vector case
            id = mycontains(Test,sprintf('%d',i)) & mycontains(Trial,sprintf('%d',j));
            if ~sum(id)  %  empty
                idvu(k) = false; k = k+1;   continue;
            end
            [~, ss{k}] = assem2d(Th,Coef(id),Test(id),Trial(id),Vhij,quadOrder);
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
if strcmpi(Test, 'v.val')
    trf = eye(nSpace); f = Coef;
    for i = 1:nSpace
        Coef = @(pz) f(pz)*trf(:, i);  Test = sprintf('v%d.val',i);
        Fi = assem2d(Th,Coef,Test,[],Vh{i},quadOrder);
        if i==1, id = 1 : NNdofv(1);  end
        if i>=2, id = NNdofv(i-1) + (1:NNdofv(i)); end
        ff(id) = ff(id) + Fi;
    end
    output = ff;  return;
end

% Test -->  'v1.val' --> {'v1.val'}
if ~iscell(Test), Coef = {Coef};  Test = {Test}; end

% Test --> {'v1.val', 'v3.val'}
for s = 1:length(Test)
    Coefv = Coef{s}; Testv = Test{s};
    for i = 1:nSpace
        str = sprintf('v%d.val',i);  % vi.val
        if strcmpi(Testv, str)
            Fi = assem2d(Th,Coefv,Testv,[],Vh{i},quadOrder);
            if i==1, id = 1 : NNdofv(1);  end
            if i>=2, id = NNdofv(i-1) + (1:NNdofv(i)); end
            ff(id) = ff(id) + Fi;
        end
    end
end
output = ff;

