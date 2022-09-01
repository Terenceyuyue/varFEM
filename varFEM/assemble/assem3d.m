function varargout = assem3d(Th,Coef,Test,Trial,Vh,quadOrder)
%%assem3d assembles bilinear or linear forms of scalar type in 3-D
%
% The key subroutine of int2d.m
% FE space:  Vh --> (v,u) --> ( Vh{1}, Vh{2} )
%

%% Preparation for the input
if nargin == 5, quadOrder = 3; end

% " v = u "
if iscell(Vh) && length(Vh)==1 % {'P1'}
    Vh = repmat(Vh, 1, 2);
end
if ~iscell(Vh) % 'P1'
    Vh = repmat( {Vh}, 1, 2);
end

if ~iscell(Coef), Coef = {Coef}; end 
if ~isempty(Test) && ~iscell(Test), Test = {Test}; end
if ~isempty(Trial) && ~iscell(Trial), Trial = {Trial}; end

% Guass-Quadrature
[lambda,weight] = quadpts3(quadOrder); 
weight = weight(:)'; % must be a row vector
ng = length(weight);
% volume
node = Th.node; elem3 = Th.elem3; NT = size(elem3,1);
z1 = node(elem3(:,1),:); 
z2 = node(elem3(:,2),:); 
z3 = node(elem3(:,3),:); 
z4 = node(elem3(:,4),:);
volume = simplexvolume(node,elem3); % volume > 0

%% Sparse assembly index
% elementwise d.o.f.s
[elem2dofv,Ndofv,NNdofv] = dof3d(Th,Vh{1}); % v
if strcmpi(Vh{1},Vh{2}) % u
    elem2dofu = elem2dofv;
    Ndofu = Ndofv;
    NNdofu = NNdofv;
else
    [elem2dofu,Ndofu,NNdofu] = dof3d(Th,Vh{2});
end
% assembly index
ii = reshape(repmat(elem2dofv, Ndofu,1), [], 1);
jj = repmat(elem2dofu(:), Ndofv, 1);

%% Coef matrix
CoefMat = cell(length(Coef),1);
for ss = 1:length(Coef)
    cf = Coef{ss};

    % case 1: cf = uh is a finite element function
    if isnumeric(cf) && length(cf)==NNdofu
        phi = Base3d('u.val',node,elem3,Vh{2},quadOrder); % phi1,phi2,phi3,phi4
        cc = zeros(NT,ng);
        for ib = 1:length(phi)
            ui = cf(elem2dofu(:,ib));
            cc = cc + repmat(ui,1,ng).*phi{ib};   % ui*phi_i
        end
    end

    % case 2: coef matrix
    if isnumeric(cf) && size(cf,1)==NT && size(cf,2)==ng
        cc = cf;
    end

    % case 3: function
    if isnumeric(cf) && length(cf)==1 % a constant
        cf = @(pz) cf+0*pz(:,1);
    end
    if ~isnumeric(cf) 
        cc = zeros(NT,ng);
        for p = 1:ng
            pz = lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3 + lambda(p,4)*z4;
            cc(:,p) = cf(pz);
        end
    end

    CoefMat{ss} = cc;
end

%% Bilinear form 
if ~isempty(Trial)
    K = zeros(NT,Ndofv*Ndofu);
    for ss = 1:length(Test)  % ss-th component of bilinear form
        % Test and Trial functions (val, dx, dy, dz, grad)
        v = Test{ss}; u = Trial{ss};
        vbase = Base3d(v,node,elem3,Vh{1},quadOrder); % v1,v2,v3,v4
        ubase = Base3d(u,node,elem3,Vh{2},quadOrder); % u1,u2,u3,u4
        
        % Coef matrix
        cc = CoefMat{ss};
        if mycontains(v,'.grad')  % of course u = u.grad
            cgrad = zeros(NT,3*ng);
            cgrad(:,1:3:end) = cc;  
            cgrad(:,2:3:end) = cc;
            cgrad(:,3:3:end) = cc; 
            cc = cgrad;
        end

        % Weight matrix
        ww = repmat(weight,NT,1);
        if mycontains(v,'.grad')  % of course u = u.grad
            ww = zeros(1,3*ng);
            ww(1:3:end) = weight; 
            ww(2:3:end) = weight;
            ww(3:3:end) = weight;
            ww = repmat(ww,NT,1);
        end

        % Stiffness matrix
        s = 1;
        for i = 1:Ndofv
            for j = 1:Ndofu
                vi = vbase{i}; uj = ubase{j};
                K(:,s) =  K(:,s) + volume.*sum(ww.*cc.*vi.*uj,2);
                s = s+1;
            end
        end
    end
    varargout{1} = sparse(ii,jj,K(:),NNdofv,NNdofu); % kk
    % K(:) is used in the assembly of vectorized FEM. This is why we assemble the
    % matrix in the form
    %      K = [k11, k12, k13, k14, k21, k22, k23, k24, ...].
    varargout{2} = K(:);  
    return;  % The remaining code will be neglected.
end

%% Linear form
if isempty(Trial) && ~isempty(Test)
    % Test functions (v.val, v.dx, v.dy, v.dz)
    F = zeros(NT,Ndofv); % straighten
    % Weight matrix
    ww = repmat(weight,NT,1);
    for ss = 1:length(Test)
        v = Test{ss};
        vbase = Base3d(v,node,elem3,Vh{1},quadOrder); % v1,v2,v3,v4

        % Coef matrix
        cc = CoefMat{ss};

        % Load vector
        for j = 1:Ndofv
            vj = vbase{j};
            F(:,j) = F(:,j) + volume.*sum(ww.*cc.*vj,2);  % [f*phi1, f*phi2, f*phi3]
        end
    end
    varargout{1} = accumarray(elem2dofv(:), F(:),[NNdofv 1]); % ff
    varargout{2} = F(:);
    return;
end

%% Integral
if isempty(Trial) && isempty(Test)
    % I = \int_\Omega Coef dx
    %   = assem3d(Th,Coef,[],[],Vh,quadOrder)
    % Weight matrix
    ww = repmat(weight,NT,1);
    I = 0;
    for ss = 1:length(Coef)
        cc = CoefMat{ss};
        I = I + dot(volume, sum(ww.*cc,2));
    end
    varargout{1} = I;
end