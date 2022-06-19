function w = Base2d_P1(wStr,node,elem,quadOrder)

% Gauss quadrature rule
[lambda,weight] = quadpts(quadOrder);  
ng = length(weight);  
NT = size(elem,1);

% gradbasis
if mycontains(wStr,'.dx') || mycontains(wStr,'.dy') || mycontains(wStr,'.grad')
    Dlambda = gradbasis(node,elem);
    Dlambda1 = Dlambda(1:NT,:,1);
    Dlambda2 = Dlambda(1:NT,:,2);
    Dlambda3 = Dlambda(1:NT,:,3);
    Dlambdax = [Dlambda1(:,1), Dlambda2(:,1), Dlambda3(:,1)];
    Dlambday = [Dlambda1(:,2), Dlambda2(:,2), Dlambda3(:,2)];
end

%% u.val
if mycontains(wStr,'.val')
    w1 = repmat(lambda(:,1)',NT,1); % phi1 at zp, p = 1,2,...
    w2 = repmat(lambda(:,2)',NT,1);
    w3 = repmat(lambda(:,3)',NT,1);
    w = {w1,w2,w3};
    return;
end
%% u.dx
if mycontains(wStr,'.dx')
    w1 = Dlambdax(:,1);  w1 = repmat(w1,1,ng);
    w2 = Dlambdax(:,2);  w2 = repmat(w2,1,ng);
    w3 = Dlambdax(:,3);  w3 = repmat(w3,1,ng);
    w = {w1,w2,w3};
    return;
end
%% u.dy
if mycontains(wStr,'.dy')
    w1 = Dlambday(:,1);  w1 = repmat(w1,1,ng);
    w2 = Dlambday(:,2);  w2 = repmat(w2,1,ng);
    w3 = Dlambday(:,3);  w3 = repmat(w3,1,ng);
    w = {w1,w2,w3};
    return;
end
%% u.grad
if mycontains(wStr,'.grad')
    w1 = Dlambda1;  w1 = repmat(w1,1,ng);
    w2 = Dlambda2;  w2 = repmat(w2,1,ng);
    w3 = Dlambda3;  w3 = repmat(w3,1,ng);
    w = {w1,w2,w3};
    return;
end