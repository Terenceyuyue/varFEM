function varargout = assem1d(Th,Coef,Test,Trial,Vh,quadOrder)
%%assem1d  assembles bilinear or linear forms of scalar type
%
% The key subroutine of int1d.m
% FE space: Vh --> (v,u) --> ( Vh{1}, Vh{2} )
%

%% Preparations
if nargin == 5, quadOrder = 3; end

% " v = u "
if iscell(Vh) && length(Vh)==1 % {'P1'}
    Vh = repmat(Vh, 1, 2);
end
if ~iscell(Vh) % 'P1'
    Vh = repmat( {Vh}, 1, 2);
end

if ~isempty(Test) && ~iscell(Test), Test = {Test}; end
if ~isempty(Trial) && ~iscell(Trial), Trial = {Trial}; end

% 1-D mesh info
node = Th.node;
if size(node,2)==1
    elem1d = Th.elem1d;
end
if size(node,2)==2 && isfield(Th,'on')
    on = Th.on;
    elem1d = Th.bdEdgeType{on};
    Th.elem1d = elem1d;
end
if isfield(Th,'elem1d'), elem1d = Th.elem1d; end  % provided by user
nel = size(elem1d,1);

% Guass-Quadrature
[lambda,weight] = quadpts1(quadOrder);
weight = weight(:)'; % must be a row vector
ng = length(weight);
% length
za = node(elem1d(:,1),:); zb = node(elem1d(:,2),:);
h = sqrt(sum((zb-za).^2,2));

%% Sparse assembly index
% elementwise d.o.f.s
[elem2dofv,Ndofv,NNdofv] = dof1d(Th,Vh{1}); % v
if strcmpi(Vh{1},Vh{2}) % u
    elem2dofu = elem2dofv;
    Ndofu = Ndofv;
    NNdofu = NNdofv;
else
    [elem2dofu,Ndofu,NNdofu] = dof1d(Th,Vh{2});
end

% assembly index
ii = reshape(repmat(elem2dofv, Ndofu,1), [], 1);
jj = repmat(elem2dofu(:), Ndofv, 1);

%% Coef matrix
if ~iscell(Coef), Coef = {Coef}; end % linear form
CoefMat = cell(length(Coef),1);
for ss = 1:length(Coef)
    cf = Coef{ss};

    % case 1: cf = uh is a finite element function ( only for dofs )
    if isnumeric(cf) && size(cf,1)==NNdofu && size(cf,2)==1
        phi = Base1d('u.val',node,elem1d,Vh{2},quadOrder); % phi1,phi2
        cc = zeros(nel,ng);
        for ib = 1:length(phi)
            ui = cf(elem2dofu(:,ib));
            cc = cc + repmat(ui,1,ng).*phi{ib};   % ui*phi_i(zp)
        end
    end

    % case 2: coef matrix
    if isnumeric(cf) && size(cf,1)==nel && size(cf,2)==ng
        cc = cf;
    end

    % case 3: function
    if isnumeric(cf) && length(cf)==1 % a constant
        cf = @(pz) cf+0*pz(:,1);
    end
    if ~isnumeric(cf)
        cc = zeros(nel,ng);
        for p = 1:ng
            pz = lambda(p,1)*za + lambda(p,2)*zb;
            cc(:,p) = cf(pz);
        end
    end

    CoefMat{ss} = cc;
end

%% Bilinear form
if ~isempty(Trial)
    K = zeros(nel,Ndofv*Ndofu); % K = [k11, k12, k21, k22]
    for ss = 1:length(Test)  % ss-th component of bilinear form

        % Test and Trial functions (.val, .dx, ...)
        v = Test{ss}; u = Trial{ss};
        vbase = Base1d(v,node,elem1d,Vh{1},quadOrder); % v1,v2
        ubase = Base1d(u,node,elem1d,Vh{2},quadOrder); % u1,u2

        % Coef matrix
        cc = CoefMat{ss};

        % Weight matrix
        ww = repmat(weight,nel,1);

        % Stiffness matrix
        s = 1;
        for i = 1:Ndofv
            for j = 1:Ndofu
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

%% Linear form
if isempty(Trial) && ~isempty(Test)
    % Test functions (v.val, v.dx, v.dy)
    F = zeros(nel,Ndofv); % straighten
    % Weight matrix
    ww = repmat(weight,nel,1);
    for ss = 1:length(Test)
        v = Test{ss};
        vbase = Base1d(v,node,elem1d,Vh{1},quadOrder); % v1,v2

        % Coef matrix
        cc = CoefMat{ss};

        % Load vector
        for j = 1:Ndofv
            vj = vbase{j};
            F(:,j) = F(:,j) + h.*sum(ww.*cc.*vj,2);  % [f*phi1, f*phi2]
        end
    end
    varargout{1} = accumarray(elem2dofv(:), F(:),[NNdofv 1]); % ff
    varargout{2} = F(:);
    return;
end

%% Integral
if isempty(Trial) && isempty(Test)
    % I = \int_\Omega Coef dx
    %   = assem1d(Th,Coef,[],[],Vh,quadOrder)
    % Weight matrix
    ww = repmat(weight,nel,1);
    I = 0;
    for ss = 1:length(Coef)
        cc = CoefMat{ss};
        I = I + dot(h, sum(ww.*cc,2));
    end
    varargout{1} = I;
end