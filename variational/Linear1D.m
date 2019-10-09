function ff = Linear1D(Th,f)

Test = {'v.val'};

node = Th.node; elem = Th.elem;
N = size(node,1); xa = node(elem(:,1)); xb = node(elem(:,2)); 

[r, w] = GaussQuad(4,0,1); % on [0,1];
xx = xa+(xb-xa)*r';  % quadrature points on all elements
J = Jacobian1D(xa,xb);  % Jacobian on all elements

v = Test{1}; [v1,v2] = Refbase1D(v);
cf = f(xx);  fe = w.*[v1(r),v2(r)]; Fref = cf*fe;
cv = Trans1D(v,xa,xb); F = cv.*J.*Fref;

ff = accumarray(elem(:), F(:),[N 1]);
