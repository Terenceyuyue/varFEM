function Cmat = getMat1d(fun,Th,quadOrder)
if nargin==2, quadOrder = 3; end

node = Th.node; 
nel = size(Th.elem1D,1);

% Guass-Quadrature
[lambda,weight] = quadpts1(quadOrder); ng = length(weight);

z1 = node(Th.elem1D(:,1),:); z2 = node(Th.elem1D(:,2),:);
Cmat = zeros(nel,ng);

% function handle is pde.Du
z0 = z1(1,:); v0 = fun(z0);
if length(v0)>1
    Du = fun;
    % nvec
    e = z1-z2;  he = sqrt(sum(e.^2,2));
    nvec = [-e(:,2)./he, e(:,1)./he];
    % Coef matrix: gN
    for p = 1:ng
        pz = lambda(p,1)*z1 + lambda(p,2)*z2;
        Cmat(:,p) = sum(Du(pz).*nvec,2);
    end
    return;  % The remaining code will be neglected.
end

% function handle is not pde.Du
for p = 1:ng
    pz = lambda(p,1)*z1 + lambda(p,2)*z2;
    Cmat(:,p) = fun(pz);
end

