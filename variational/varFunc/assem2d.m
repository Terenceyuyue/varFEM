function varargout = assem2d(Th,Coef,Test,Trial,Vh,quadOrder)
%% ASSEM2D assembles bilinear or linear forms of scalar type
%
% The key subroutine of int2d.m
% FE space:  Vh --> (v,u) --> ( Vh{1}, Vh{2} )
% Note: For programming in terms of block matrix, use assem2d.m rather than
% int2d.m
%

%% Preparation for the input
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

% 2D mesh information
node = Th.node;
if size(node,2) == 2
    elem = Th.elem;
else
    elem = Th.elem2D;   % add in programming ( 3-D boundary )
end
NT = size(elem,1);

% Guass-Quadrature
[lambda,weight] = quadpts(quadOrder); 
weight = weight(:)'; % must be a row vector
nG = length(weight);
% area
z1 = node(elem(:,1),:); z2 = node(elem(:,2),:); z3 = node(elem(:,3),:);
area = simplexvolume(node,elem);

%% Sparse assembling index of Pk-Lagrange element 
% elementwise d.o.f.s
[elem2dofv,Ndofv,NNdofv] = dof2d(Th,Vh{1}); % v
[elem2dofu,Ndofu,NNdofu] = dof2d(Th,Vh{2}); % u

% assembling index
ii = reshape(repmat(elem2dofv, Ndofu,1), [], 1);
jj = repmat(elem2dofu(:), Ndofv, 1);

%% Bilinear form 
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
            cgrad = zeros(NT,2*nG);
            cgrad(:,1:2:end) = cc;  
            cgrad(:,2:2:end) = cc; 
            cc = cgrad;
        end
        % Weight matrix
        ww = repmat(weight,NT,1);
        if mycontains(v,'.grad')  % of course u = u.grad
            ww = zeros(1,2*nG);
            ww(1:2:end) = weight; 
            ww(2:2:end) = weight;
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
    % K(:) is used in the assembly of vectorized FEM. This is why we assemble the
    % matrix in the form
    %      K = [k11, k12, k13, k21, k22, k23, k31, k32, k33].
    varargout{2} = K(:);  
    return;  % The remaining code will be neglected.
end

%% Linear form 
% % isempty(Trial) 
% Test functions (v.val)
v = Test;
vbase = Base2D(v,node,elem,Vh{1},quadOrder); % v1,v2,v3
% Coef matrix
if isnumeric(Coef)&&size(Coef,2)==nG, cc = Coef; end % matrix of size NT*nG
if length(Coef)==1   % function handle or a constant
    cf = Coef;
    if isnumeric(cf), cf = @(p) cf+0*p(:,1); end % a constant
    cc = zeros(NT,nG);
    for p = 1:nG
        pz = lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
        cc(:,p) = cf(pz);
    end
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

