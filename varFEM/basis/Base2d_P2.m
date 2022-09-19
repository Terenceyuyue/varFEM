function w = Base2d_P2(wStr,node,elem,quadOrder)

% Gauss quadrature rule
[lambda,weight] = quadpts(quadOrder);  
ng = length(weight);  
NT = size(elem,1);

dot = strfind(wStr,'.');
wStr = wStr(dot:end);

% gradbasis
if ~strcmpi(wStr,'.val')
    Dlambda = gradbasis(node,elem);
    Dlambda1 = Dlambda(1:NT,:,1);
    Dlambda2 = Dlambda(1:NT,:,2);
    Dlambda3 = Dlambda(1:NT,:,3);
    Dlambdax = [Dlambda1(:,1), Dlambda2(:,1), Dlambda3(:,1)];
    Dlambday = [Dlambda1(:,2), Dlambda2(:,2), Dlambda3(:,2)];
end


%% u.val
if strcmpi(wStr,'.val')
    w1 = lambda(:,1)'.*(2*lambda(:,1)'-1);
    w2 = lambda(:,2)'.*(2*lambda(:,2)'-1);
    w3 = lambda(:,3)'.*(2*lambda(:,3)'-1);
    w4 = 4*lambda(:,2)'.*lambda(:,3)';
    w5 = 4*lambda(:,1)'.*lambda(:,3)';
    w6 = 4*lambda(:,1)'.*lambda(:,2)';
    w1 = repmat(w1,NT,1); w2 = repmat(w2,NT,1); w3 = repmat(w3,NT,1);
    w4 = repmat(w4,NT,1); w5 = repmat(w5,NT,1); w6 = repmat(w6,NT,1);
    w = {w1,w2,w3,w4,w5,w6};
    return;
end
%% u.dx
if strcmpi(wStr,'.dx')
    [w1,w2,w3,w4,w5,w6] = deal(zeros(NT,ng));
    for p = 1:ng
        w1(:,p) = Dlambdax(:,1)*(4*lambda(p,1)-1);
        w2(:,p) = Dlambdax(:,2)*(4*lambda(p,2)-1);
        w3(:,p) = Dlambdax(:,3)*(4*lambda(p,3)-1);
        w4(:,p) = 4*(Dlambdax(:,2)*lambda(p,3) + Dlambdax(:,3)*lambda(p,2));
        w5(:,p) = 4*(Dlambdax(:,1)*lambda(p,3) + Dlambdax(:,3)*lambda(p,1));
        w6(:,p) = 4*(Dlambdax(:,1)*lambda(p,2) + Dlambdax(:,2)*lambda(p,1));
    end
    w = {w1,w2,w3,w4,w5,w6};
    return;
end
%% u.dy
if strcmpi(wStr,'.dy')
    [w1,w2,w3,w4,w5,w6] = deal(zeros(NT,ng));
    for p = 1:ng
        w1(:,p) = Dlambday(:,1)*(4*lambda(p,1)-1);
        w2(:,p) = Dlambday(:,2)*(4*lambda(p,2)-1);
        w3(:,p) = Dlambday(:,3)*(4*lambda(p,3)-1);
        w4(:,p) = 4*(Dlambday(:,2)*lambda(p,3) + Dlambday(:,3)*lambda(p,2));
        w5(:,p) = 4*(Dlambday(:,1)*lambda(p,3) + Dlambday(:,3)*lambda(p,1));
        w6(:,p) = 4*(Dlambday(:,1)*lambda(p,2) + Dlambday(:,2)*lambda(p,1));
    end
    w = {w1,w2,w3,w4,w5,w6};
    return;
end
%% u.grad
if strcmpi(wStr,'.grad')
    [w1,w2,w3,w4,w5,w6] = deal(zeros(NT,2*ng));
    for p = 1:ng
        w1(:,2*p-1:2*p) = Dlambda1*(4*lambda(p,1)-1);
        w2(:,2*p-1:2*p) = Dlambda2*(4*lambda(p,2)-1);
        w3(:,2*p-1:2*p) = Dlambda3*(4*lambda(p,3)-1);
        w4(:,2*p-1:2*p) = 4*(Dlambda2*lambda(p,3) + Dlambda3*lambda(p,2));
        w5(:,2*p-1:2*p) = 4*(Dlambda1*lambda(p,3) + Dlambda3*lambda(p,1));
        w6(:,2*p-1:2*p) = 4*(Dlambda1*lambda(p,2) + Dlambda2*lambda(p,1));
    end
    w = {w1,w2,w3,w4,w5,w6};
    return;
end
%% u.dxx
if strcmpi(wStr,'.dxx')
    w1 = repmat(4*Dlambdax(:,1).^2,1,ng);
    w2 = repmat(4*Dlambdax(:,2).^2,1,ng);
    w3 = repmat(4*Dlambdax(:,3).^2,1,ng);
    w4 = repmat(8*(Dlambdax(:,2).*Dlambdax(:,3)),1,ng);
    w5 = repmat(8*(Dlambdax(:,1).*Dlambdax(:,3)),1,ng);
    w6 = repmat(8*(Dlambdax(:,1).*Dlambdax(:,2)),1,ng);
    w = {w1,w2,w3,w4,w5,w6};
    return;
end
%% u.dyy
if strcmpi(wStr,'.dyy')
    w1 = repmat(4*Dlambday(:,1).^2,1,ng);
    w2 = repmat(4*Dlambday(:,2).^2,1,ng);
    w3 = repmat(4*Dlambday(:,3).^2,1,ng);
    w4 = repmat(8*(Dlambday(:,2).*Dlambday(:,3)),1,ng);
    w5 = repmat(8*(Dlambday(:,1).*Dlambday(:,3)),1,ng);
    w6 = repmat(8*(Dlambday(:,1).*Dlambday(:,2)),1,ng);
    w = {w1,w2,w3,w4,w5,w6};
    return;
end
%% u.dyy
if strcmpi(wStr,'.dxy') || strcmpi(wStr,'.dyx')
    w1 = repmat(4*Dlambdax(:,1).*Dlambday(:,1),1,ng);
    w2 = repmat(4*Dlambdax(:,2).*Dlambday(:,2),1,ng);
    w3 = repmat(4*Dlambdax(:,3).*Dlambday(:,3),1,ng);
    w4 = repmat(4*(Dlambdax(:,2).*Dlambday(:,3) + Dlambday(:,2).*Dlambdax(:,3)),1,ng);
    w5 = repmat(4*(Dlambdax(:,1).*Dlambday(:,3) + Dlambday(:,1).*Dlambdax(:,3)),1,ng);
    w6 = repmat(4*(Dlambdax(:,1).*Dlambday(:,2) + Dlambday(:,1).*Dlambdax(:,2)),1,ng);
    w = {w1,w2,w3,w4,w5,w6};
    return;
end