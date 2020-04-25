function varargout = assem1d(Th,Coef,Test,Trial,Vh,quadOrder)

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

% 1D mesh information
node = Th.node;
elem1D = Th.elem1D;  nel = size(elem1D,1);

% Guass-Quadrature
[lambda,weight] = quadpts1(quadOrder); ng = length(weight);
% length
za = node(elem1D(:,1),:); zb = node(elem1D(:,2),:);
h = sqrt(sum((zb-za).^2,2));

%% ------------------- Sparse assembling index ---------------------
% elementwise d.o.f.s
[elem2dofv,ndofv,NNdofv] = dof1d(Th,Vh{1}); % v
[elem2dofu,ndofu,NNdofu] = dof1d(Th,Vh{2}); % u

% assembling index
nnz = nel*ndofv*ndofu;
ii = zeros(nnz,1); jj = zeros(nnz,1); id = 0;
for i = 1:ndofv
    for j = 1:ndofu
        ii(id+1:id+nel) = elem2dofv(:,i);   
        jj(id+1:id+nel) = elem2dofu(:,j);   
        id = id + nel;
    end
end

%% ------------------------- Bilinear form -------------------------
if ~isempty(Trial)
    K = zeros(nel,ndofv*ndofu);
    for ss = 1:length(Test)  % ss-th component of bilinear form
        % Test and Trial functions (val, dx, ...)
        v = Test{ss}; u = Trial{ss};
        vbase = Base1D(v,node,elem1D,Vh{1},quadOrder); % v1,v2
        ubase = Base1D(u,node,elem1D,Vh{2},quadOrder); % u1,u2
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
        for i = 1:ndofv
            for j = 1:ndofu
                vi = vbase{i}; uj = ubase{j};
                K(:,s) =  K(:,s) + h.*sum(ww.*cc.*vi.*uj,2);
                s = s+1;
            end
        end
    end
    varargout{1} = sparse(ii,jj,K(:),NNdofv,NNdofu); % kk
    varargout{2} = K(:); 
    return;  % The remaining code will be neglected.
end

%% ------------------------- Linear form -------------------------
% % if isempty(Trial)
% Test functions (v.val)
v = Test;
vbase = Base1D(v,node,elem1D,Vh{1},quadOrder); % v1,v2
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
F = zeros(nel,ndofv); % straighten
for j = 1:ndofv
    vj = vbase{j};
    F(:,j) = h.*sum(ww.*cc.*vj,2);  % [f*phi1, f*phi2]
end
varargout{1} = accumarray(elem2dofv(:), F(:),[NNdofv 1]); % ff
varargout{2} = F(:);