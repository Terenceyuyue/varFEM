function output = int1d(Th,Coef,Trial,Test,feSpace,quadOrder)

if nargin == 4, feSpace = 'P1'; quadOrder = 3; end % default: P1
if nargin == 5, quadOrder = 3; end

% mesh information
node = Th.node;
elem1D = Th.elem1D;  nel = size(elem1D,1);

% Guass-Quadrature
[lambda,weight] = quadpts1(quadOrder); ng = length(weight);
% length
za = node(elem1D(:,1),:); zb = node(elem1D(:,2),:);
h = sqrt(sum((zb-za).^2,2));

%% ------------------- Sparse assembling index ---------------------
% elementwise d.o.f.s
[elem2dof,ndof,NNdof] = dof1d(Th,feSpace);

% assembling index
nnz = nel*ndof^2;
ii = zeros(nnz,1); jj = zeros(nnz,1); id = 0;
for i = 1:ndof
    for j = 1:ndof
        ii(id+1:id+nel) = elem2dof(:,i);   % zi
        jj(id+1:id+nel) = elem2dof(:,j);   % zj
        id = id + nel;
    end
end

%% ------------------------- Bilinear form -------------------------
if ~isempty(Trial)
    K = zeros(nel,ndof^2);
    for ss = 1:length(Test)  % ss-th component of bilinear form
        % Test and Trial functions (val, dx)
        v = Test{ss}; u = Trial{ss};
        vbase = Base1D(v,node,elem1D,feSpace,quadOrder); % v1,v2
        ubase = Base1D(u,node,elem1D,feSpace,quadOrder); % u1,u2
        % Coef matrix
        cf = Coef{ss};
        if isnumeric(cf) && length(cf)==1, cf = @(pz) cf+0*pz(:,1); end  % a constant, x = pz(:,1)
        cc = zeros(nel,ng);
        for p = 1:ng
            pz = lambda(p,1)*za + lambda(p,2)*zb;
            cc(:,p) = cf(pz);
        end
        % Weight matrix
        ww = repmat(weight,nel,1);
        % Stiffness matrix
        s = 1;
        for i = 1:ndof
            for j = 1:ndof
                vi = vbase{i}; uj = ubase{j};
                K(:,s) =  K(:,s) + h.*sum(ww.*cc.*vi.*uj,2);
                s = s+1;
            end
        end
    end
    output = sparse(ii,jj,K(:),NNdof,NNdof); % kk
    return;  % The remaining code will be neglected.
end

%% ------------------------- Linear form -------------------------
% % if isempty(Trial)
% Test functions (v.val)
v = Test;
vbase = Base1D(v,node,elem1D,feSpace,quadOrder); % v1,v2
% Coef matrix
if isnumeric(Coef)&&size(Coef,2)==ng, cc = Coef; end % matrix of size NT*ng
if length(Coef)==1  % function handle or a constant
    cf = Coef;
    if isnumeric(cf), cf = @(p) cf + 0*p(:,1); end
    cc = zeros(nel,ng);
    for p = 1:ng
        pz = lambda(p,1)*za + lambda(p,2)*zb;
        cc(:,p) = cf(pz);
    end
end
% Weight matrix
ww = repmat(weight,nel,1);
% Load vector
F = zeros(nel,ndof); % straighten
for j = 1:ndof
    vj = vbase{j};
    F(:,j) = h.*sum(ww.*cc.*vj,2);  % [f*phi1, f*phi2]
end
output = accumarray(elem2dof(:), F(:),[NNdof 1]); % ff


