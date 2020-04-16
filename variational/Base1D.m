function w = Base1D(wStr,node,elem1D,feSpace,quadOrder)

if nargin == 3, feSpace = []; quadOrder = 3; end % default: P1
if nargin == 4, quadOrder = 3; end

wStr = lower(wStr); % lowercase string
nel = size(elem1D,1);

% length
za = node(elem1D(:,1),:); zb = node(elem1D(:,2),:);
h = sqrt(sum((zb-za).^2,2));

% Gauss quadrature rule
[lambda,weight] = quadpts1(quadOrder); ng = length(weight);

% derivatives of bases
Dlambda1 = -1./h;
Dlambda2 = 1./h;

%% P1-Lagrange
if  strcmpi(feSpace, 'P1') || isempty(feSpace)
    % u.val
    if mycontains(wStr,'.val') 
        w1 = repmat(lambda(:,1)',nel,1); % phi1 at xp, p = 1,2,...
        w2 = repmat(lambda(:,2)',nel,1);
    end
    % u.dx
    if mycontains(wStr,'.dx') 
        w1 = Dlambda1;  w1 = repmat(w1,1,ng);
        w2 = Dlambda2;  w2 = repmat(w2,1,ng);
    end
    w = {w1,w2};
end

%% P2-Lagrange
if  strcmpi(feSpace, 'P2')
    % u.val
    if mycontains(wStr,'.val')
        w1 = lambda(:,1)'.*(2*lambda(:,1)'-1);  w1 = repmat(w1,nel,1);
        w2 = lambda(:,2)'.*(2*lambda(:,2)'-1);  w2 = repmat(w2,nel,1);
        w3 = 4*lambda(:,1)'.*lambda(:,2)';      w3 = repmat(w3,nel,1);
    end
    % u.dx
    if mycontains(wStr,'.dx')
        w1 = zeros(nel,ng);  w2 = w1; w3 = w1;
        for p = 1:ng
            w1(:,p) = (4*lambda(p,1)-1)*Dlambda1;
            w2(:,p) = (4*lambda(p,2)-1)*Dlambda2;
            w3(:,p) = 4*(lambda(p,1)*Dlambda2+lambda(p,2)*Dlambda1);
        end
    end
    w = {w1,w2,w3};
end

%% P3-Lagrange
if  strcmpi(feSpace, 'P3')
    % u.val
    if mycontains(wStr,'.val')
        w1 = 0.5*(3*lambda(:,1)-1).*(3*lambda(:,1)-2).*lambda(:,1);
        w2 = 0.5*(3*lambda(:,2)-1).*(3*lambda(:,2)-2).*lambda(:,2);
        w3 = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,1)-1);
        w4 = 9/2*lambda(:,2).*lambda(:,1).*(3*lambda(:,2)-1);
        w1 = repmat(w1',nel,1);      w2 = repmat(w2',nel,1);
        w3 = repmat(w3',nel,1);      w4 = repmat(w4',nel,1);
    end
    % u.dx
    if mycontains(wStr,'.dx')
        w1 = zeros(nel,ng);  w2 = w1; w3 = w1;
        for p = 1:ng
            w1(:,p) = (27/2*lambda(p,1)^2-9*lambda(p,1)+1).*Dlambda1;
            w2(:,p) = (27/2*lambda(p,2)^2-9*lambda(p,2)+1).*Dlambda2;
            w3(:,p) = 9/2*(lambda(p,1)*(3*lambda(p,1)-1).*Dlambda2 ...
                    + lambda(p,2)*(6*lambda(p,1)-1).*Dlambda1);
            w4(:,p) = 9/2*(lambda(p,2)*(3*lambda(p,2)-1).*Dlambda1 ...
                    + lambda(p,1)*(6*lambda(p,2)-1).*Dlambda2);
        end
    end
    w = {w1,w2,w3,w4};
end
