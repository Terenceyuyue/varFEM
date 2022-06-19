function Cmat = interpEdgeMat(fun,Th,quadOrder)
%%interpEdgeMat returns the matrix for numerical integration on boundary edges
%
%  This function is designed to deal with terms containing normal
%  vector, for example, fun = pde.Du*n
%

if nargin==2, quadOrder = 3; end

node = Th.node;
if isfield(Th,'on')
    on = Th.on;
    elem1d = Th.bdEdgeType{on};
end
if isfield(Th,'elem1d'), elem1d = Th.elem1d; end
nel = size(elem1d,1); % 1-D mesh on the boundary

% Guass-Quadrature
[lambda,weight] = quadpts1(quadOrder); ng = length(weight);

za = node(elem1d(:,1),:);
zb = node(elem1d(:,2),:);
Cmat = zeros(nel,ng);  % size: nel * ng

% function handle is pde.Du
z0 = za(1,:); 
v0 = fun(z0);
if length(v0)>1
    Du = fun;
    % nvec
    e = za-zb;  
    he = sqrt(sum(e.^2,2));
    nvec = [-e(:,2)./he, e(:,1)./he];
    % Coef matrix: gN
    for p = 1:ng
        pz = lambda(p,1)*za + lambda(p,2)*zb;
        Cmat(:,p) = sum(Du(pz).*nvec,2);
    end
    return;  % The remaining code will be neglected.
end

% function handle is not pde.Du
for p = 1:ng
    pz = lambda(p,1)*za + lambda(p,2)*zb;
    Cmat(:,p) = fun(pz);
end