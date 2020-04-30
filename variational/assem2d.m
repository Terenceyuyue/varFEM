function varargout = assem2d(Th,Coef,Test,Trial,Vh,quadOrder)

% Vh --> (v,u) --> ( Vh{1}, Vh{2} )

%% ---------------- Preparation for the input -----------------
% default para
if nargin==3 % Trial = [] --> linear form
    Trial = []; Vh = {'P1','P1'}; quadOrder = 3; 
end
if nargin == 4  % default: (v,u) --> (P1,P1)
    Vh = {'P1','P1'}; quadOrder = 3; 
end 
if nargin == 5, quadOrder = 3; end

% " v = u "
if iscell(Vh) && length(Vh)==1 % {'P1'}
    Vh = repmat(Vh, 1, 2); 
end 
if ~iscell(Vh) % 'P1'
    Vh = repmat( {Vh}, 1, 2); 
end  

% Guass-Quadrature
[lambda,weight] = quadpts(quadOrder); nG = length(weight);
% area
node = Th.node; elem = Th.elem; NT = size(elem,1);
z1 = node(elem(:,1),:); z2 = node(elem(:,2),:); z3 = node(elem(:,3),:);
xi = [z2(:,1)-z3(:,1), z3(:,1)-z1(:,1), z1(:,1)-z2(:,1)];
eta = [z2(:,2)-z3(:,2), z3(:,2)-z1(:,2), z1(:,2)-z2(:,2)];
area = 0.5*(xi(:,1).*eta(:,2)-xi(:,2).*eta(:,1));

%% ------ Sparse assembling index of Pk-Lagrange element ------
% elementwise d.o.f.s
[elem2dofv,Ndofv,NNdofv] = dof2d(Th,Vh{1}); % v
[elem2dofu,Ndofu,NNdofu] = dof2d(Th,Vh{2}); % u

% assembling index
nnz = NT*Ndofv*Ndofu;
ii = zeros(nnz,1); jj = zeros(nnz,1); id = 0;
for i = 1:Ndofv
    for j = 1:Ndofu
        ii(id+1:id+NT) = elem2dofv(:,i);   
        jj(id+1:id+NT) = elem2dofu(:,j);   
        id = id + NT;
    end
end

%% ------------------------- Bilinear form -------------------------
if ~isempty(Trial)
    K = zeros(NT,Ndofv*Ndofu);
    for ss = 1:length(Test)  % ss-th component of bilinear form
        % Test and Trial functions (val, dx, dy, grad)
        v = Test{ss}; u = Trial{ss};
        vbase = Base2D(v,node,elem,Vh{1},quadOrder); % v1,v2,v3
        ubase = Base2D(u,node,elem,Vh{2},quadOrder); % u1,u2,u3
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
        for i = 1:Ndofv
            for j = 1:Ndofu
                vi = vbase{i}; uj = ubase{j};
                K(:,s) =  K(:,s) + area.*sum(ww.*cc.*vi.*uj,2);
                s = s+1;
            end
        end
    end
    varargout{1} = sparse(ii,jj,K(:),NNdofv,NNdofu); % kk
    varargout{2} = K(:);
    return;  % The remaining code will be neglected.
end

%% ------------------------- Linear form -------------------------
% % isempty(Trial) 
% Test functions (v.val)
v = Test;
vbase = Base2D(v,node,elem,Vh{1},quadOrder); % v1,v2,v3
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
F = zeros(NT,Ndofv); % straighten
for j = 1:Ndofv
    vj = vbase{j};
    F(:,j) = area.*sum(ww.*cc.*vj,2);  % [f*phi1, f*phi2, f*phi3]
end
varargout{1} = accumarray(elem2dofv(:), F(:),[NNdofv 1]); % ff
varargout{2} = F(:);

