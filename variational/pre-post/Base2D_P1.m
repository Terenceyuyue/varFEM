% This is not a separate matlab script storing the base matrix of P1 Lagrange element.
% see Base2D.m
%
% Copyright (C) Terence Yu.

%% u.val
if mycontains(wStr,'.val')
    w1 = repmat(lambda(:,1)',NT,1); % phi1 at zp, p = 1,2,...
    w2 = repmat(lambda(:,2)',NT,1);
    w3 = repmat(lambda(:,3)',NT,1);
end
%% u.dx
if mycontains(wStr,'.dx')
    w1 = Dlambdax(:,1);  w1 = repmat(w1,1,nQuad);
    w2 = Dlambdax(:,2);  w2 = repmat(w2,1,nQuad);
    w3 = Dlambdax(:,3);  w3 = repmat(w3,1,nQuad);
end
%% u.dy
if mycontains(wStr,'.dy')
    w1 = Dlambday(:,1);  w1 = repmat(w1,1,nQuad);
    w2 = Dlambday(:,2);  w2 = repmat(w2,1,nQuad);
    w3 = Dlambday(:,3);  w3 = repmat(w3,1,nQuad);
end
%% u.grad
if mycontains(wStr,'.grad')
    w1 = Dlambda1;  w1 = repmat(w1,1,nQuad);
    w2 = Dlambda2;  w2 = repmat(w2,1,nQuad);
    w3 = Dlambda3;  w3 = repmat(w3,1,nQuad);
end

w = {w1,w2,w3};