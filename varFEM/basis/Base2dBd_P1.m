function w = Base2dBd_P1(wStr,node,elem,quadOrder)

% Gauss quadrature rule for elementwise edges
[lambdaBd,weightBd] = quadptsBd(quadOrder);  
nG = length(weightBd);  % nG = 3*ng
NT = size(elem,1);

dot = strfind(wStr,'.');
wStr = wStr(dot:end);

% gradbasis
if ~strcmpi(wStr,'.val')
    Dlambda = gradbasis(node,elem);  % constant on K
    Dlambda1 = Dlambda(1:NT,:,1);
    Dlambda2 = Dlambda(1:NT,:,2);
    Dlambda3 = Dlambda(1:NT,:,3);
    Dlambdax = [Dlambda1(:,1), Dlambda2(:,1), Dlambda3(:,1)];
    Dlambday = [Dlambda1(:,2), Dlambda2(:,2), Dlambda3(:,2)];
end

%% u.val
if strcmpi(wStr,'.val')
    w1 = repmat(lambdaBd(:,1)',NT,1); % phi1 at zp, p = 1,2,...
    w2 = repmat(lambdaBd(:,2)',NT,1);
    w3 = repmat(lambdaBd(:,3)',NT,1);
    w = {w1,w2,w3};
    return;
end
%% u.dx
if strcmpi(wStr,'.dx')
    w1 = Dlambdax(:,1);  w1 = repmat(w1,1,nG);
    w2 = Dlambdax(:,2);  w2 = repmat(w2,1,nG);
    w3 = Dlambdax(:,3);  w3 = repmat(w3,1,nG);
    w = {w1,w2,w3};
    return;
end
%% u.dy
if strcmpi(wStr,'.dy')
    w1 = Dlambday(:,1);  w1 = repmat(w1,1,nG);
    w2 = Dlambday(:,2);  w2 = repmat(w2,1,nG);
    w3 = Dlambday(:,3);  w3 = repmat(w3,1,nG);
    w = {w1,w2,w3};
    return;
end
%% u.grad
if strcmpi(wStr,'.grad')
    w1 = Dlambda1;  w1 = repmat(w1,1,nG);
    w2 = Dlambda2;  w2 = repmat(w2,1,nG);
    w3 = Dlambda3;  w3 = repmat(w3,1,nG);
    w = {w1,w2,w3};
    return;
end