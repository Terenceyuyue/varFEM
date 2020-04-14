function varargout = int2d(Th,Coef,Trial,Test,feSpace,quadOrder)

if nargin == 4, feSpace = 'P1'; quadOrder = 3; end % default: P1
if nargin == 5, quadOrder = 3; end

node = Th.node; elem = Th.elem;
NT = size(elem,1);

% Guass-Quadrature
[lambda,weight] = quadpts(quadOrder); nG = length(weight);
% area
z1 = node(elem(:,1),:); z2 = node(elem(:,2),:); z3 = node(elem(:,3),:);
xi = [z2(:,1)-z3(:,1), z3(:,1)-z1(:,1), z1(:,1)-z2(:,1)];
eta = [z2(:,2)-z3(:,2), z3(:,2)-z1(:,2), z1(:,2)-z2(:,2)];
area = 0.5*(xi(:,1).*eta(:,2)-xi(:,2).*eta(:,1));

%% ------ Sparse assembling index of Pk-Lagrange element ------
% elementwise d.o.f.s
[elem2dof,Ndof,NNdof] = dof2d(Th,feSpace);

% assembling index
nnz = NT*Ndof^2;
ii = zeros(nnz,1); jj = zeros(nnz,1); id = 0;
for i = 1:Ndof
    for j = 1:Ndof
        ii(id+1:id+NT) = elem2dof(:,i);   % zi
        jj(id+1:id+NT) = elem2dof(:,j);   % zj
        id = id + NT;
    end
end

%% ------------------------- Bilinear form -------------------------
if ~isempty(Trial)
    K = zeros(NT,Ndof^2);
    for ss = 1:length(Test)  % ss-th component of bilinear form
        % Test and Trial functions (val, dx, dy, grad)
        v = Test{ss}; u = Trial{ss};
        vbase = Base2D(v,node,elem,feSpace,quadOrder); % v1,v2,v3
        ubase = Base2D(u,node,elem,feSpace,quadOrder); % u1,u2,u3
        % Coef matrix
        cf = Coef{ss};
        if isnumeric(cf) && length(cf)==1, cf = @(pz) cf+0*pz(:,1); end
        cc = zeros(NT,nG);
        for p = 1:nG
            pz = lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
            cc(:,p) = cf(pz);
        end
        if mycontains(v,'.grad')  % of course u = u.grad
            cgrad = ones(NT,2*nG);
            cgrad(:,1:2:end) = cc;  cgrad(:,2:2:end) = cc; cc = cgrad;
        end
        % Weight matrix
        ww = repmat(weight,NT,1);
        if mycontains(v,'.grad')  % of course u = u.grad
            ww = zeros(1,2*nG);
            ww(1:2:end) = weight; ww(2:2:end) = weight;
            ww = repmat(ww,NT,1);
        end
        % Stiffness matrix
        s = 1;
        for i = 1:Ndof
            for j = 1:Ndof
                vi = vbase{i}; uj = ubase{j};
                K(:,s) =  K(:,s) + area.*sum(ww.*cc.*vi.*uj,2);
                s = s+1;
            end
        end
    end
    varargout{1} = sparse(ii,jj,K(:),NNdof,NNdof); % kk
    varargout{2} = K(:);
    return;  % The remaining code will be neglected.
end

%% ------------------------- Linear form -------------------------
% % isempty(Trial) 
% Test functions (v.val)
v = Test;
vbase = Base2D(v,node,elem,feSpace,quadOrder); % v1,v2,v3
% Coef matrix
f = Coef;
if isnumeric(f) && length(f)==1, f = @(pz) f+0*pz(:,1); end % a constant
cc = zeros(NT,nG);
for p = 1:nG
    pz = lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
    cc(:,p) = f(pz);
end
% Weight matrix
ww = repmat(weight,NT,1);
% Load vector
F = zeros(NT,Ndof); % straighten
for j = 1:Ndof
    vj = vbase{j};
    F(:,j) = area.*sum(ww.*cc.*vj,2);  % [f*phi1, f*phi2, f*phi3]
end
varargout{1} = accumarray(elem2dof(:), F(:),[NNdof 1]); % ff
varargout{2} = F(:);

