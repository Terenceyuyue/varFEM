function output = int2dvec(Th,Coef,Trial,Test,feSpace,quadOrder)

if nargin == 4, feSpace = 'P1'; quadOrder = 3; end % default: P1
if nargin == 5, quadOrder = 3; end

%% ------------------- extended [Coef,Trial,Test] ------------------
% % example
% Coef = { 1, 1, 0.5 }; % may be function handles
% Trial = {'u1.dx', 'u2.dy', 'u1.dy + u2.dx'};
% Test  = {'v1.dx', 'v2.dy', 'v1.dy + v2.dx'};

if ~isempty(Trial)
    cc = cell(1,length(Coef)); uu = cc; vv = cc;
    for ss = 1:length(Test)
        stru = Trial{ss}; strv = Test{ss};
        [cs, us, vs] = getvarForm(stru, strv);
        cf = Coef{ss};
        for i = 1:length(cs)
            if isnumeric(cf), cs{i} = cf*cs{i}; end
            if isa(cf, 'function_handle'), cs{i} = @(p) cf(p)*cs{i}; end
        end
        cc{ss} = cs; uu{ss} = us; vv{ss} = vs;
    end
    
    Coef  = horzcat(cc{:});
    Trial = horzcat(uu{:});
    Test  = horzcat(vv{:});
end
%% --------- Sparse assembling index of Pk-Lagrange element --------
% elementwise d.o.f.s
[elem2dof,Ndof,NNdof] = dof2d(Th,feSpace);

% assembling index
NT = size(Th.elem,1);
nnz = NT*Ndof^2;
ii = zeros(nnz,1); jj = zeros(nnz,1); id = 0;
for i = 1:Ndof
    for j = 1:Ndof
        ii(id+1:id+NT) = elem2dof(:,i);   % zi
        jj(id+1:id+NT) = elem2dof(:,j);   % zj
        id = id + NT;
    end
end

ii11 = ii;       jj11 = jj;  ii12 = ii;       jj12 = jj+NNdof;
ii21 = ii+NNdof; jj21 = jj;  ii22 = ii+NNdof; jj22 = jj+NNdof;

ii = [ii11; ii12; ii21; ii22];
jj = [jj11; jj12; jj21; jj22];

%% ------------------------- Bilinear form -------------------------
if ~isempty(Trial)
    ss = zeros(nnz,4); k = 1;
    % (vi,uj), i,j = 1,2
    for i = 1:2
        for j = 1:2
            id = mycontains(Test,sprintf('%d',i)) & mycontains(Trial,sprintf('%d',j));
            if ~sum(id), k = k+1;   continue;  end   %  empty 
            [~, ss(:,k)] = int2d(Th,Coef(id),Trial(id),Test(id),feSpace,quadOrder);
            k = k+1;
        end
    end
    ss = ss(:);
    output = sparse(ii,jj,ss,2*NNdof,2*NNdof);
    return;  % The remaining code will be neglected.
end

%% -------------------------- Linear form --------------------------
F = zeros(2*NNdof, 1);
if strcmpi(Test, 'v1.val') % v1.val
    F1 = int2d(Th,Coef,[],Test,feSpace,quadOrder);
    F(1:NNdof) = F(1:NNdof) + F1;
end
if strcmpi(Test, 'v2.val') % v2.val
    F2 = int2d(Th,Coef,[],Test,feSpace,quadOrder);
    F(NNdof+1:end) = F(NNdof+1:end) + F2;
end
if strcmpi(Test, 'v.val')  % v.val = [v1.val, v2.val]
    trf = eye(2); f = Coef;
    Coef = @(pz) f(pz)*trf(:, 1);  Test = 'v1.val';
    F1 = int2d(Th,Coef,[],Test,feSpace,quadOrder);
    Coef = @(pz) f(pz)*trf(:, 2);  Test = 'v2.val';
    F2 = int2d(Th,Coef,[],Test,feSpace,quadOrder);
    F = [F1; F2];
end
output = F;



