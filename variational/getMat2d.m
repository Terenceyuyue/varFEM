function Cmat = getMat2d(fun,Th,quadOrder)
if nargin==2, quadOrder = 3; end

node = Th.node; elem = Th.elem; NT = size(elem,1);
z1 = node(elem(:,1),:); z2 = node(elem(:,2),:); z3 = node(elem(:,3),:);

% Guass-Quadrature
[lambda,weight] = quadpts(quadOrder); nG = length(weight);


Cmat = zeros(NT,nG);  % size: NT * nG

for p = 1:nG
    pz = lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
    Cmat(:,p) = fun(pz);
end

