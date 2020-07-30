function w = Base2D(wStr,node,elem,Vh,quadOrder)
%% BASE2D returns base matrix for numerical integration
%
%  w = {w1,w2, ..., wn}, where n is the number of the basis functions
%  wi is a matrix of size NT*nG, where NT and nG are the numbers of
%  2D elements and quadrature points


if nargin == 3, Vh = 'P1'; quadOrder = 3; end % default: P1
if nargin == 4, quadOrder = 3; end

%% Quadrature and gradbasis
wStr = lower(wStr); % lowercase string
NT = size(elem,1);

% Gauss quadrature rule
[lambda,weight] = quadpts(quadOrder); nQuad = length(weight);

% gradbasis
Dlambda = gradbasis(node,elem);
Dlambda1 = Dlambda(1:NT,:,1);
Dlambda2 = Dlambda(1:NT,:,2);
Dlambda3 = Dlambda(1:NT,:,3);
Dlambdax = [Dlambda1(:,1), Dlambda2(:,1), Dlambda3(:,1)];
Dlambday = [Dlambda1(:,2), Dlambda2(:,2), Dlambda3(:,2)];

%% P1-Lagrange
if  strcmpi(Vh, 'P1')
    % u.val
    if mycontains(wStr,'.val')
        w1 = repmat(lambda(:,1)',NT,1); % phi1 at zp, p = 1,2,...
        w2 = repmat(lambda(:,2)',NT,1);
        w3 = repmat(lambda(:,3)',NT,1);
    end
    % u.dx
    if mycontains(wStr,'.dx')
        w1 = Dlambdax(:,1);  w1 = repmat(w1,1,nQuad);
        w2 = Dlambdax(:,2);  w2 = repmat(w2,1,nQuad);
        w3 = Dlambdax(:,3);  w3 = repmat(w3,1,nQuad);
    end
    % u.dy
    if mycontains(wStr,'.dy')
        w1 = Dlambday(:,1);  w1 = repmat(w1,1,nQuad);
        w2 = Dlambday(:,2);  w2 = repmat(w2,1,nQuad);
        w3 = Dlambday(:,3);  w3 = repmat(w3,1,nQuad);
    end
    % u.grad
    if mycontains(wStr,'.grad')
        w1 = Dlambda1;  w1 = repmat(w1,1,nQuad);
        w2 = Dlambda2;  w2 = repmat(w2,1,nQuad);
        w3 = Dlambda3;  w3 = repmat(w3,1,nQuad);
    end
    
    w = {w1,w2,w3};
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    % u.val
    if mycontains(wStr,'.val')
        w1 = lambda(:,1)'.*(2*lambda(:,1)'-1);
        w2 = lambda(:,2)'.*(2*lambda(:,2)'-1);
        w3 = lambda(:,3)'.*(2*lambda(:,3)'-1);
        w4 = 4*lambda(:,2)'.*lambda(:,3)';
        w5 = 4*lambda(:,1)'.*lambda(:,3)';
        w6 = 4*lambda(:,1)'.*lambda(:,2)';
        w1 = repmat(w1,NT,1); w2 = repmat(w2,NT,1); w3 = repmat(w3,NT,1);
        w4 = repmat(w4,NT,1); w5 = repmat(w5,NT,1); w6 = repmat(w6,NT,1);
    end
    % u.dx
    if mycontains(wStr,'.dx')
        w1 = zeros(NT,nQuad); w2 = w1; w3 = w1; w4 = w1; w5 = w1; w6 = w1;
        for p = 1:nQuad
            w1(:,p) = Dlambdax(:,1)*(4*lambda(p,1)-1);
            w2(:,p) = Dlambdax(:,2)*(4*lambda(p,2)-1);
            w3(:,p) = Dlambdax(:,3)*(4*lambda(p,3)-1);
            w4(:,p) = 4*(Dlambdax(:,2)*lambda(p,3) + Dlambdax(:,3)*lambda(p,2));
            w5(:,p) = 4*(Dlambdax(:,1)*lambda(p,3) + Dlambdax(:,3)*lambda(p,1));
            w6(:,p) = 4*(Dlambdax(:,1)*lambda(p,2) + Dlambdax(:,2)*lambda(p,1));
        end
    end
    % u.dy
    if mycontains(wStr,'.dy')
        w1 = zeros(NT,nQuad); w2 = w1; w3 = w1; w4 = w1; w5 = w1; w6 = w1;
        for p = 1:nQuad
            w1(:,p) = Dlambday(:,1)*(4*lambda(p,1)-1);
            w2(:,p) = Dlambday(:,2)*(4*lambda(p,2)-1);
            w3(:,p) = Dlambday(:,3)*(4*lambda(p,3)-1);
            w4(:,p) = 4*(Dlambday(:,2)*lambda(p,3) + Dlambday(:,3)*lambda(p,2));
            w5(:,p) = 4*(Dlambday(:,1)*lambda(p,3) + Dlambday(:,3)*lambda(p,1));
            w6(:,p) = 4*(Dlambday(:,1)*lambda(p,2) + Dlambday(:,2)*lambda(p,1));
        end
    end
    % u.grad
    if mycontains(wStr,'.grad')
        w1 = zeros(NT,2*nQuad); w2 = w1; w3 = w1; w4 = w1; w5 = w1; w6 = w1;
        for p = 1:nQuad
            w1(:,2*p-1:2*p) = Dlambda1*(4*lambda(p,1)-1);
            w2(:,2*p-1:2*p) = Dlambda2*(4*lambda(p,2)-1);
            w3(:,2*p-1:2*p) = Dlambda3*(4*lambda(p,3)-1);
            w4(:,2*p-1:2*p) = 4*(Dlambda2*lambda(p,3) + Dlambda3*lambda(p,2));
            w5(:,2*p-1:2*p) = 4*(Dlambda1*lambda(p,3) + Dlambda3*lambda(p,1));
            w6(:,2*p-1:2*p) = 4*(Dlambda1*lambda(p,2) + Dlambda2*lambda(p,1));
        end
    end
    
    w = {w1,w2,w3,w4,w5,w6};
end

%% P3-Lagrange
%     % 3                                 
%     % |  \                                
%     % 5    7                              
%     % |     \                    
%     % 8  10  4                            
%     % |        \                          
%     % 1- 6 - 9 - 2                       
% The order of the d.o.f.s is shown on the left.

if  strcmpi(Vh, 'P3')
    % u.val
    if mycontains(wStr,'.val')
        w1 = 0.5*(3*lambda(:,1)-1).*(3*lambda(:,1)-2).*lambda(:,1);
        w2 = 0.5*(3*lambda(:,2)-1).*(3*lambda(:,2)-2).*lambda(:,2);
        w3 = 0.5*(3*lambda(:,3)-1).*(3*lambda(:,3)-2).*lambda(:,3);
        w4 = 9/2*lambda(:,3).*lambda(:,2).*(3*lambda(:,2)-1);
        w5 = 9/2*lambda(:,1).*lambda(:,3).*(3*lambda(:,3)-1);
        w6 = 9/2*lambda(:,2).*lambda(:,1).*(3*lambda(:,1)-1);
        w7 = 9/2*lambda(:,2).*lambda(:,3).*(3*lambda(:,3)-1);
        w8 = 9/2*lambda(:,3).*lambda(:,1).*(3*lambda(:,1)-1);
        w9 = 9/2*lambda(:,2).*lambda(:,1).*(3*lambda(:,2)-1);
        w10 = 27*lambda(:,1).*lambda(:,2).*lambda(:,3);
        w1 = repmat(w1',NT,1); w2 = repmat(w2',NT,1); w3 = repmat(w3',NT,1);
        w4 = repmat(w4',NT,1); w5 = repmat(w5',NT,1); w6 = repmat(w6',NT,1);
        w7 = repmat(w7',NT,1); w8 = repmat(w8',NT,1); w9 = repmat(w9',NT,1);
        w10 = repmat(w10',NT,1);
    end
    % u.dx
    if mycontains(wStr,'.dx')
        w1 = zeros(NT,nQuad); w2 = w1; w3 = w1; w4 = w1; w5 = w1; w6 = w1; w7 = w1; w8 = w1; w9 = w1; w10 = w1;
        for p = 1:nQuad
            w1(:,p) = 0.5*(3*Dlambdax(:,1)).*(3*lambda(p,1)-2).*lambda(p,1) ...
                + 0.5*(3*lambda(p,1)-1).*(3*Dlambdax(:,1)).*lambda(p,1) ...
                + 0.5*(3*lambda(p,1)-1).*(3*lambda(p,1)-2).*Dlambdax(:,1);
            w2(:,p) = 0.5*(3*Dlambdax(:,2)).*(3*lambda(p,2)-2).*lambda(p,2) ...
                + 0.5*(3*lambda(p,2)-1).*(3*Dlambdax(:,2)).*lambda(p,2) ...
                + 0.5*(3*lambda(p,2)-1).*(3*lambda(p,2)-2).*Dlambdax(:,2);
            w3(:,p) = 0.5*(3*Dlambdax(:,3)).*(3*lambda(p,3)-2).*lambda(p,3) ...
                + 0.5*(3*lambda(p,3)-1).*(3*Dlambdax(:,3)).*lambda(p,3) ...
                + 0.5*(3*lambda(p,3)-1).*(3*lambda(p,3)-2).*Dlambdax(:,3);
            w4(:,p) = 9/2*Dlambdax(:,3).*lambda(p,2).*(3*lambda(p,2)-1) ...
                + 9/2*lambda(p,3).*Dlambdax(:,2).*(3*lambda(p,2)-1) ...
                + 9/2*lambda(p,3).*lambda(p,2).*(3*Dlambdax(:,2));
            w5(:,p) = 9/2*Dlambdax(:,1).*lambda(p,3).*(3*lambda(p,3)-1) ...
                + 9/2*lambda(p,1).*Dlambdax(:,3).*(3*lambda(p,3)-1) ...
                + 9/2*lambda(p,1).*lambda(p,3).*(3*Dlambdax(:,3));
            w6(:,p) = 9/2*Dlambdax(:,1).*lambda(p,2).*(3*lambda(p,1)-1) ...
                + 9/2*lambda(p,1).*Dlambdax(:,2).*(3*lambda(p,1)-1) ...
                + 9/2*lambda(p,1).*lambda(p,2).*(3*Dlambdax(:,1));
            w7(:,p) = 9/2*Dlambdax(:,3).*lambda(p,2).*(3*lambda(p,3)-1) ...
                + 9/2*lambda(p,3).*Dlambdax(:,2).*(3*lambda(p,3)-1) ...
                + 9/2*lambda(p,3).*lambda(p,2).*(3*Dlambdax(:,3));
            w8(:,p) = 9/2*Dlambdax(:,1).*lambda(p,3).*(3*lambda(p,1)-1) ...
                + 9/2*lambda(p,1).*Dlambdax(:,3).*(3*lambda(p,1)-1) ...
                + 9/2*lambda(p,1).*lambda(p,3).*(3*Dlambdax(:,1));
            w9(:,p) = 9/2*Dlambdax(:,1).*lambda(p,2).*(3*lambda(p,2)-1) ...
                + 9/2*lambda(p,1).*Dlambdax(:,2).*(3*lambda(p,2)-1) ...
                + 9/2*lambda(p,1).*lambda(p,2).*(3*Dlambdax(:,2));
            w10(:,p) = 27*Dlambdax(:,1).*lambda(p,2).*lambda(p,3) ...
                + 27*lambda(p,1).*Dlambdax(:,2).*lambda(p,3) ...
                + 27*lambda(p,1).*lambda(p,2).*Dlambdax(:,3);
        end
    end
    % u.dy
    if mycontains(wStr,'.dy')
        w1 = zeros(NT,nQuad); w2 = w1; w3 = w1; w4 = w1; w5 = w1; w6 = w1;
        for p = 1:nQuad
            w1(:,p) = 0.5*(3*Dlambday(:,1)).*(3*lambda(p,1)-2).*lambda(p,1) ...
                + 0.5*(3*lambda(p,1)-1).*(3*Dlambday(:,1)).*lambda(p,1) ...
                + 0.5*(3*lambda(p,1)-1).*(3*lambda(p,1)-2).*Dlambday(:,1);
            w2(:,p) = 0.5*(3*Dlambday(:,2)).*(3*lambda(p,2)-2).*lambda(p,2) ...
                + 0.5*(3*lambda(p,2)-1).*(3*Dlambday(:,2)).*lambda(p,2) ...
                + 0.5*(3*lambda(p,2)-1).*(3*lambda(p,2)-2).*Dlambday(:,2);
            w3(:,p) = 0.5*(3*Dlambday(:,3)).*(3*lambda(p,3)-2).*lambda(p,3) ...
                + 0.5*(3*lambda(p,3)-1).*(3*Dlambday(:,3)).*lambda(p,3) ...
                + 0.5*(3*lambda(p,3)-1).*(3*lambda(p,3)-2).*Dlambday(:,3);
            w4(:,p) = 9/2*Dlambday(:,3).*lambda(p,2).*(3*lambda(p,2)-1) ...
                + 9/2*lambda(p,3).*Dlambday(:,2).*(3*lambda(p,2)-1) ...
                + 9/2*lambda(p,3).*lambda(p,2).*(3*Dlambday(:,2));
            w5(:,p) = 9/2*Dlambday(:,1).*lambda(p,3).*(3*lambda(p,3)-1) ...
                + 9/2*lambda(p,1).*Dlambday(:,3).*(3*lambda(p,3)-1) ...
                + 9/2*lambda(p,1).*lambda(p,3).*(3*Dlambday(:,3));
            w6(:,p) = 9/2*Dlambday(:,1).*lambda(p,2).*(3*lambda(p,1)-1) ...
                + 9/2*lambda(p,1).*Dlambday(:,2).*(3*lambda(p,1)-1) ...
                + 9/2*lambda(p,1).*lambda(p,2).*(3*Dlambday(:,1));
            w7(:,p) = 9/2*Dlambday(:,3).*lambda(p,2).*(3*lambda(p,3)-1) ...
                + 9/2*lambda(p,3).*Dlambday(:,2).*(3*lambda(p,3)-1) ...
                + 9/2*lambda(p,3).*lambda(p,2).*(3*Dlambday(:,3));
            w8(:,p) = 9/2*Dlambday(:,1).*lambda(p,3).*(3*lambda(p,1)-1) ...
                + 9/2*lambda(p,1).*Dlambday(:,3).*(3*lambda(p,1)-1) ...
                + 9/2*lambda(p,1).*lambda(p,3).*(3*Dlambday(:,1));
            w9(:,p) = 9/2*Dlambday(:,1).*lambda(p,2).*(3*lambda(p,2)-1) ...
                + 9/2*lambda(p,1).*Dlambday(:,2).*(3*lambda(p,2)-1) ...
                + 9/2*lambda(p,1).*lambda(p,2).*(3*Dlambday(:,2));
            w10(:,p) = 27*Dlambday(:,1).*lambda(p,2).*lambda(p,3) ...
                + 27*lambda(p,1).*Dlambday(:,2).*lambda(p,3) ...
                + 27*lambda(p,1).*lambda(p,2).*Dlambday(:,3);
        end
    end
    % u.grad
    if mycontains(wStr,'.grad')
        w1 = zeros(NT,2*nQuad); w2 = w1; w3 = w1; w4 = w1; w5 = w1; w6 = w1; w7 = w1; w8 = w1; w9 = w1; w10 = w1;
        for p = 1:nQuad
            w1(:,2*p-1:2*p) = 0.5*(3*Dlambda1).*(3*lambda(p,1)-2).*lambda(p,1) ...
                + 0.5*(3*lambda(p,1)-1).*(3*Dlambda1).*lambda(p,1) ...
                + 0.5*(3*lambda(p,1)-1).*(3*lambda(p,1)-2).*Dlambda1;
            w2(:,2*p-1:2*p) = 0.5*(3*Dlambda2).*(3*lambda(p,2)-2).*lambda(p,2) ...
                + 0.5*(3*lambda(p,2)-1).*(3*Dlambda2).*lambda(p,2) ...
                + 0.5*(3*lambda(p,2)-1).*(3*lambda(p,2)-2).*Dlambda2;
            w3(:,2*p-1:2*p) = 0.5*(3*Dlambda3).*(3*lambda(p,3)-2).*lambda(p,3) ...
                + 0.5*(3*lambda(p,3)-1).*(3*Dlambda3).*lambda(p,3) ...
                + 0.5*(3*lambda(p,3)-1).*(3*lambda(p,3)-2).*Dlambda3;
            w4(:,2*p-1:2*p) = 9/2*Dlambda3.*lambda(p,2).*(3*lambda(p,2)-1) ...
                + 9/2*lambda(p,3).*Dlambda2.*(3*lambda(p,2)-1) ...
                + 9/2*lambda(p,3).*lambda(p,2).*(3*Dlambda2);
            w5(:,2*p-1:2*p) = 9/2*Dlambda1.*lambda(p,3).*(3*lambda(p,3)-1) ...
                + 9/2*lambda(p,1).*Dlambda3.*(3*lambda(p,3)-1) ...
                + 9/2*lambda(p,1).*lambda(p,3).*(3*Dlambda3);
            w6(:,2*p-1:2*p) = 9/2*Dlambda1.*lambda(p,2).*(3*lambda(p,1)-1) ...
                + 9/2*lambda(p,1).*Dlambda2.*(3*lambda(p,1)-1) ...
                + 9/2*lambda(p,1).*lambda(p,2).*(3*Dlambda1);
            w7(:,2*p-1:2*p) = 9/2*Dlambda3.*lambda(p,2).*(3*lambda(p,3)-1) ...
                + 9/2*lambda(p,3).*Dlambda2.*(3*lambda(p,3)-1) ...
                + 9/2*lambda(p,3).*lambda(p,2).*(3*Dlambda3);
            w8(:,2*p-1:2*p) = 9/2*Dlambda1.*lambda(p,3).*(3*lambda(p,1)-1) ...
                + 9/2*lambda(p,1).*Dlambda3.*(3*lambda(p,1)-1) ...
                + 9/2*lambda(p,1).*lambda(p,3).*(3*Dlambda1);
            w9(:,2*p-1:2*p) = 9/2*Dlambda1.*lambda(p,2).*(3*lambda(p,2)-1) ...
                + 9/2*lambda(p,1).*Dlambda2.*(3*lambda(p,2)-1) ...
                + 9/2*lambda(p,1).*lambda(p,2).*(3*Dlambda2);
            w10(:,2*p-1:2*p) = 27*Dlambda1.*lambda(p,2).*lambda(p,3) ...
                + 27*lambda(p,1).*Dlambda2.*lambda(p,3) ...
                + 27*lambda(p,1).*lambda(p,2).*Dlambda3;
        end
    end
    
    w = {w1,w2,w3,w4,w5,w6,w7,w8,w9,w10};
end