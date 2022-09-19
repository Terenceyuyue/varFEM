function w = Base3d_P2(wStr,node,elem3,quadOrder)

% Gauss quadrature rule
[lambda,weight] = quadpts3(quadOrder);
ng = length(weight);
NT = size(elem3,1);

% gradbasis
if ~strcmpi(wStr,'.val')
    Dlambda = gradbasis3(node,elem3);
    Dlambda1 = Dlambda(1:NT,:,1);
    Dlambda2 = Dlambda(1:NT,:,2);
    Dlambda3 = Dlambda(1:NT,:,3);
    Dlambda4 = Dlambda(1:NT,:,4);
    Dlambdax = Dlambda(1:NT,1,1:4);
    Dlambday = Dlambda(1:NT,2,1:4);
    Dlambdaz = Dlambda(1:NT,3,1:4);
end

%% u.val
if strcmpi(wStr,'.val')
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
    w = {w1,w2,w3,w4,w5,w6,w7,w8,w9,w10};
    return;
end
%% u.dx
if strcmpi(wStr,'.dx')
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10] = deal(zeros(NT,ng));
    for p = 1:ng
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
    w = {w1,w2,w3,w4,w5,w6,w7,w8,w9,w10};
    return;
end
%% u.dy
if strcmpi(wStr,'.dy')
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10] = deal(zeros(NT,ng));
    for p = 1:ng
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
    w = {w1,w2,w3,w4,w5,w6,w7,w8,w9,w10};
    return;
end

%% u.dz
if strcmpi(wStr,'.dz')
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10] = deal(zeros(NT,ng));
    for p = 1:ng
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
    w = {w1,w2,w3,w4,w5,w6,w7,w8,w9,w10};
    return;
end

%% u.grad
if strcmpi(wStr,'.grad')
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10] = deal(zeros(NT,3*ng));
    for p = 1:ng
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
    w = {w1,w2,w3,w4,w5,w6,w7,w8,w9,w10};
    return;
end