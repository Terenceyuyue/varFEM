function w = Base3D(wStr,node,elem,Vh,quadOrder)
%% BASE3D returns base matrix for numerical integration
%
%  w = {w1,w2, ..., wn}, where n is the number of the basis functions
%  wi is a matrix of size NT*nG, where NT and nG are the numbers of
%  3D elements and quadrature points
%
% Copyright (C) Terence Yu.

if nargin == 3, Vh = 'P1'; quadOrder = 3; end % default: P1
if nargin == 4, quadOrder = 3; end

%% Quadrature and gradbasis
wStr = lower(wStr); % lowercase string
NT = size(elem,1);

% Gauss quadrature rule
[lambda,weight] = quadpts3(quadOrder); nQuad = length(weight);

% gradbasis
Dlambda = gradbasis3(node,elem);
Dlambda1 = Dlambda(1:NT,:,1);
Dlambda2 = Dlambda(1:NT,:,2);
Dlambda3 = Dlambda(1:NT,:,3);
Dlambda4 = Dlambda(1:NT,:,4);
Dlambdax = Dlambda(1:NT,1,1:4);
Dlambday = Dlambda(1:NT,2,1:4);
Dlambdaz = Dlambda(1:NT,3,1:4);

%% P1-Lagrange
if  strcmpi(Vh, 'P1')
    % u.val
    if mycontains(wStr,'.val')
        w1 = repmat(lambda(:,1)',NT,1); % phi1 at zp, p = 1,2,...
        w2 = repmat(lambda(:,2)',NT,1);
        w3 = repmat(lambda(:,3)',NT,1);
        w4 = repmat(lambda(:,4)',NT,1);
    end
    % u.dx
    if mycontains(wStr,'.dx')
        w1 = Dlambdax(:,1);  w1 = repmat(w1,1,nQuad);
        w2 = Dlambdax(:,2);  w2 = repmat(w2,1,nQuad);
        w3 = Dlambdax(:,3);  w3 = repmat(w3,1,nQuad);
        w4 = Dlambdax(:,4);  w4 = repmat(w4,1,nQuad);
    end
    % u.dy
    if mycontains(wStr,'.dy')
        w1 = Dlambday(:,1);  w1 = repmat(w1,1,nQuad);
        w2 = Dlambday(:,2);  w2 = repmat(w2,1,nQuad);
        w3 = Dlambday(:,3);  w3 = repmat(w3,1,nQuad);
        w4 = Dlambday(:,4);  w4 = repmat(w4,1,nQuad);
    end
    % u.dz
    if mycontains(wStr,'.dz')
        w1 = Dlambdaz(:,1);  w1 = repmat(w1,1,nQuad);
        w2 = Dlambdaz(:,2);  w2 = repmat(w2,1,nQuad);
        w3 = Dlambdaz(:,3);  w3 = repmat(w3,1,nQuad);
        w4 = Dlambdaz(:,4);  w4 = repmat(w4,1,nQuad);
    end
    % u.grad
    if mycontains(wStr,'.grad')
        w1 = Dlambda1;  w1 = repmat(w1,1,nQuad);
        w2 = Dlambda2;  w2 = repmat(w2,1,nQuad);
        w3 = Dlambda3;  w3 = repmat(w3,1,nQuad);
        w4 = Dlambda4;  w4 = repmat(w4,1,nQuad);
    end
    
    w = {w1,w2,w3,w4};
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    % u.val
    if mycontains(wStr,'.val')        
        w1 = lambda(:,1)'.*(2*lambda(:,1)'-1);  % a1
        w2 = lambda(:,2)'.*(2*lambda(:,2)'-1);  % a2
        w3 = lambda(:,3)'.*(2*lambda(:,3)'-1);  % a3
        w4 = lambda(:,4)'.*(2*lambda(:,4)'-1);  % a4
        w5 = 4*lambda(:,1)'.*lambda(:,2)';      % a12
        w6 = 4*lambda(:,1)'.*lambda(:,3)';      % a13
        w7 = 4*lambda(:,1)'.*lambda(:,4)';      % a14
        w8 = 4*lambda(:,2)'.*lambda(:,3)';      % a23
        w9 = 4*lambda(:,2)'.*lambda(:,4)';      % a23
        w10= 4*lambda(:,3)'.*lambda(:,4)';      % a34
        w1 = repmat(w1,NT,1); w2 = repmat(w2,NT,1); w3 = repmat(w3,NT,1);
        w4 = repmat(w4,NT,1); w5 = repmat(w5,NT,1); w6 = repmat(w6,NT,1);
        w7 = repmat(w7,NT,1); w8 = repmat(w8,NT,1); w9 = repmat(w9,NT,1);
        w10 = repmat(w10,NT,1);
    end
    % u.dx
    if mycontains(wStr,'.dx')
        w1 = zeros(NT,nQuad); w2 = w1; w3 = w1; w4 = w1; 
        w5 = w1; w6 = w1; w7 = w1; w8 = w1; w9 = w1; w10 = w1; 
        for p = 1:nQuad
            w1(:,p) = (4*lambda(p,1)-1).*Dlambdax(:,1);
            w2(:,p) = (4*lambda(p,2)-1).*Dlambdax(:,2);
            w3(:,p) = (4*lambda(p,3)-1).*Dlambdax(:,3);
            w4(:,p) = (4*lambda(p,4)-1).*Dlambdax(:,4);
            w5(:,p) = 4*(lambda(p,1)*Dlambdax(:,2)+lambda(p,2)*Dlambdax(:,1));
            w6(:,p) = 4*(lambda(p,1)*Dlambdax(:,3)+lambda(p,3)*Dlambdax(:,1));
            w7(:,p) = 4*(lambda(p,1)*Dlambdax(:,4)+lambda(p,4)*Dlambdax(:,1));
            w8(:,p) = 4*(lambda(p,2)*Dlambdax(:,3)+lambda(p,3)*Dlambdax(:,2));
            w9(:,p) = 4*(lambda(p,2)*Dlambdax(:,4)+lambda(p,4)*Dlambdax(:,2));
            w10(:,p) = 4*(lambda(p,3)*Dlambdax(:,4)+lambda(p,4)*Dlambdax(:,3));
        end
    end
    % u.dy
    if mycontains(wStr,'.dy')
        w1 = zeros(NT,nQuad); w2 = w1; w3 = w1; w4 = w1; 
        w5 = w1; w6 = w1; w7 = w1; w8 = w1; w9 = w1; w10 = w1; 
        for p = 1:nQuad
            w1(:,p) = (4*lambda(p,1)-1).*Dlambday(:,1);
            w2(:,p) = (4*lambda(p,2)-1).*Dlambday(:,2);
            w3(:,p) = (4*lambda(p,3)-1).*Dlambday(:,3);
            w4(:,p) = (4*lambda(p,4)-1).*Dlambday(:,4);
            w5(:,p) = 4*(lambda(p,1)*Dlambday(:,2)+lambda(p,2)*Dlambday(:,1));
            w6(:,p) = 4*(lambda(p,1)*Dlambday(:,3)+lambda(p,3)*Dlambday(:,1));
            w7(:,p) = 4*(lambda(p,1)*Dlambday(:,4)+lambda(p,4)*Dlambday(:,1));
            w8(:,p) = 4*(lambda(p,2)*Dlambday(:,3)+lambda(p,3)*Dlambday(:,2));
            w9(:,p) = 4*(lambda(p,2)*Dlambday(:,4)+lambda(p,4)*Dlambday(:,2));
            w10(:,p) = 4*(lambda(p,3)*Dlambday(:,4)+lambda(p,4)*Dlambday(:,3));
        end
    end
    % u.dz
    if mycontains(wStr,'.dz')
        w1 = zeros(NT,nQuad); w2 = w1; w3 = w1; w4 = w1; 
        w5 = w1; w6 = w1; w7 = w1; w8 = w1; w9 = w1; w10 = w1; 
        for p = 1:nQuad
            w1(:,p) = (4*lambda(p,1)-1).*Dlambdaz(:,1);
            w2(:,p) = (4*lambda(p,2)-1).*Dlambdaz(:,2);
            w3(:,p) = (4*lambda(p,3)-1).*Dlambdaz(:,3);
            w4(:,p) = (4*lambda(p,4)-1).*Dlambdaz(:,4);
            w5(:,p) = 4*(lambda(p,1)*Dlambdaz(:,2)+lambda(p,2)*Dlambdaz(:,1));
            w6(:,p) = 4*(lambda(p,1)*Dlambdaz(:,3)+lambda(p,3)*Dlambdaz(:,1));
            w7(:,p) = 4*(lambda(p,1)*Dlambdaz(:,4)+lambda(p,4)*Dlambdaz(:,1));
            w8(:,p) = 4*(lambda(p,2)*Dlambdaz(:,3)+lambda(p,3)*Dlambdaz(:,2));
            w9(:,p) = 4*(lambda(p,2)*Dlambdaz(:,4)+lambda(p,4)*Dlambdaz(:,2));
            w10(:,p) = 4*(lambda(p,3)*Dlambdaz(:,4)+lambda(p,4)*Dlambdaz(:,3));
        end
    end
    % u.grad
    if mycontains(wStr,'.grad')
        w1 = zeros(NT,3*nQuad); w2 = w1; w3 = w1; w4 = w1; w5 = w1; w6 = w1;
        w7 = w1; w8 = w1; w9 = w1; w10 = w1; 
        for p = 1:nQuad
            w1(:,3*p-2:3*p) = (4*lambda(p,1)-1).*Dlambda1;
            w2(:,3*p-2:3*p) = (4*lambda(p,2)-1).*Dlambda2;
            w3(:,3*p-2:3*p) = (4*lambda(p,3)-1).*Dlambda3;
            w4(:,3*p-2:3*p) = (4*lambda(p,4)-1).*Dlambda4;
            w5(:,3*p-2:3*p) = 4*(lambda(p,1)*Dlambda2+lambda(p,2)*Dlambda1);
            w6(:,3*p-2:3*p) = 4*(lambda(p,1)*Dlambda3+lambda(p,3)*Dlambda1);
            w7(:,3*p-2:3*p) = 4*(lambda(p,1)*Dlambda4+lambda(p,4)*Dlambda1);
            w8(:,3*p-2:3*p) = 4*(lambda(p,2)*Dlambda3+lambda(p,3)*Dlambda2);
            w9(:,3*p-2:3*p) = 4*(lambda(p,2)*Dlambda4+lambda(p,4)*Dlambda2);
            w10(:,3*p-2:3*p) = 4*(lambda(p,3)*Dlambda4+lambda(p,4)*Dlambda3);
        end
    end
    
    w = {w1,w2,w3,w4,w5,w6,w7,w8,w9,w10};
end




