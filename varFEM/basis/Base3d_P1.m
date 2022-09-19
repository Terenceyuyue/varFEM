function w = Base3d_P1(wStr,node,elem3,quadOrder)

% Gauss quadrature rule
[lambda,weight] = quadpts3(quadOrder);
ng = length(weight);
NT = size(elem3,1);

dot = strfind(wStr,'.');
wStr = wStr(dot:end);

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
    w1 = repmat(lambda(:,1)',NT,1); % phi1 at zp, p = 1,2,...
    w2 = repmat(lambda(:,2)',NT,1);
    w3 = repmat(lambda(:,3)',NT,1);
    w4 = repmat(lambda(:,4)',NT,1);
    w = {w1,w2,w3,w4};
    return;
end
%% u.dx
if strcmpi(wStr,'.dx')
    w1 = Dlambdax(:,1);  w1 = repmat(w1,1,ng);
    w2 = Dlambdax(:,2);  w2 = repmat(w2,1,ng);
    w3 = Dlambdax(:,3);  w3 = repmat(w3,1,ng);
    w4 = Dlambdax(:,4);  w4 = repmat(w4,1,ng);
    w = {w1,w2,w3,w4};
    return;
end
%% u.dy
if strcmpi(wStr,'.dy')
    w1 = Dlambday(:,1);  w1 = repmat(w1,1,ng);
    w2 = Dlambday(:,2);  w2 = repmat(w2,1,ng);
    w3 = Dlambday(:,3);  w3 = repmat(w3,1,ng);
    w4 = Dlambday(:,4);  w4 = repmat(w4,1,ng);
    w = {w1,w2,w3,w4};
    return;
end

%% u.dz
if strcmpi(wStr,'.dz')
    w1 = Dlambdaz(:,1);  w1 = repmat(w1,1,ng);
    w2 = Dlambdaz(:,2);  w2 = repmat(w2,1,ng);
    w3 = Dlambdaz(:,3);  w3 = repmat(w3,1,ng);
    w4 = Dlambdaz(:,4);  w4 = repmat(w4,1,ng);
    w = {w1,w2,w3,w4};
    return;
end
%% u.grad
if strcmpi(wStr,'.grad')
    w1 = Dlambda1;  w1 = repmat(w1,1,ng);
    w2 = Dlambda2;  w2 = repmat(w2,1,ng);
    w3 = Dlambda3;  w3 = repmat(w3,1,ng);
    w4 = Dlambda4;  w4 = repmat(w4,1,ng);
    w = {w1,w2,w3,w4};
    return;
end