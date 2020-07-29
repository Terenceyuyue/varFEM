function err = getH1error1D(node,elem1D,uh,Du,feSpace,quadOrder)

if nargin == 4, feSpace = []; quadOrder = 3; end
if nargin == 5, quadOrder = 3; end

N = size(node,1); nel = size(elem1D,1);
% Gauss quadrature rule
[lambda,weight] = quadpts1(quadOrder); ng = length(weight);
% length
z1 = node(elem1D(:,1),:); z2 = node(elem1D(:,2),:); 
h = sqrt(sum((z2-z1).^2,2));
% derivatives of bases
Dlambda1 = -1./h; 
Dlambda2 = 1./h;

%% P1-Lagrange
if strcmpi(feSpace,'P1') || isempty(feSpace)
    % elementwise d.o.f.s
    elem2dof = elem1D;
    % numerical gradient  (at p-th quadrature point)
    Duh =  uh(elem2dof(:,1)).*Dlambda1 ...
        +  uh(elem2dof(:,2)).*Dlambda2;
    % elementwise error
    err = zeros(nel,1);
    for p = 1:ng
        pz =  lambda(p,1)*z1 + lambda(p,2)*z2;
        err = err + weight(p)*sum((Du(pz)-Duh).^2,2);
    end
end

%% P2-Lagrange
if strcmpi(feSpace,'P2')
    % elementwise d.o.f.s
    elem2dof = [elem1D, (1:nel)'+N];
    % numerical gradient (at p-th quadrature point)
    err = zeros(nel,1);
    for p = 1:ng
        Dphip1 = (4*lambda(p,1)-1)*Dlambda1;
        Dphip2 = (4*lambda(p,2)-1)*Dlambda2;
        Dphip3 = 4*(lambda(p,1)*Dlambda2+lambda(p,2)*Dlambda1);
        Duh =  uh(elem2dof(:,1)).*Dphip1 ...
            +  uh(elem2dof(:,2)).*Dphip2 ...
            +  uh(elem2dof(:,3)).*Dphip3;
        % elementwise error
        pz =  lambda(p,1)*z1 + lambda(p,2)*z2;
        err = err + weight(p)*sum((Du(pz)-Duh).^2,2);
    end    
end

%% P3-Lagrange
if strcmpi(feSpace,'P3')
    % elementwise d.o.f.s
    elem2dof = [elem1D, (1:nel)'+N, (1:nel)'+N+nel];
    % numerical gradient (at p-th quadrature point)
    err = zeros(nel,1);
    for p = 1:ng
        Dphip1 = (27/2*lambda(p,1)^2-9*lambda(p,1)+1).*Dlambda1;           
        Dphip2 = (27/2*lambda(p,2)^2-9*lambda(p,2)+1).*Dlambda2;        
        Dphip3 = 9/2*(lambda(p,1)*(3*lambda(p,1)-1).*Dlambda2 ...
                    + lambda(p,2)*(6*lambda(p,1)-1).*Dlambda1);       
        Dphip4 = 9/2*(lambda(p,2)*(3*lambda(p,2)-1).*Dlambda1 ...
                    + lambda(p,1)*(6*lambda(p,2)-1).*Dlambda2);        
        Duh =  uh(elem2dof(:,1)).*Dphip1 ...
            +  uh(elem2dof(:,2)).*Dphip2 ...
            +  uh(elem2dof(:,3)).*Dphip3 ...
            +  uh(elem2dof(:,4)).*Dphip4;
        % elementwise error
        pz =  lambda(p,1)*z1 + lambda(p,2)*z2;
        err = err + weight(p)*sum((Du(pz)-Duh).^2,2);
    end    
end

%% Modification
err = h.*err;
err(isnan(err)) = 0; % singular values are excluded
err = sqrt(abs(sum(err)));