function w = Base1d(wStr,node,elem1d,Vh,quadOrder)
%% BASE1D  returns base matrix for numerical integration
%
%  w = {w1,w2, ..., wn}, where n is the number of the basis functions
%  wi is a matrix of size nel*ng, where nel and ng are the numbers of
%  1D elements and quadrature points
%

if nargin == 3, Vh = []; quadOrder = 3; end % default: P1
if nargin == 4, quadOrder = 3; end

wStr = lower(wStr); % lowercase string

%% Quadrature and gradbasis
nel = size(elem1d,1);
% length
h = sqrt(sum((node(elem1d(:,2),:)-node(elem1d(:,1),:)).^2,2));
% Gauss quadrature rule
[lambda,weight] = quadpts1(quadOrder); ng = length(weight);
% derivatives of basis functions
Dlambda1 = -1./h;
Dlambda2 = 1./h;

%% P1-Lagrange
if  strcmpi(Vh, 'P1') || isempty(Vh)
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
    return;
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
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
    return;
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3')
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
    return;
end