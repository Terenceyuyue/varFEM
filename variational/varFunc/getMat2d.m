function Cmat = getMat2d(fun,Th,quadOrder)
%% GETMAT1D returns the matrix for numerical integration in int1d
%
%  This function is designed to deal with the terms containing normal
%  vector, for example, fun = pde.Du*n
%
% Copyright (C) Terence Yu.

if nargin==2, quadOrder = 3; end

node = Th.node; 
nf = size(Th.elem2D,1); % 2-D mesh on th boundary

% Guass-Quadrature
[lambda,weight] = quadpts(quadOrder); nG = length(weight);

z1 = node(Th.elem2D(:,1),:); 
z2 = node(Th.elem2D(:,2),:); 
z3 = node(Th.elem2D(:,3),:);
Cmat = zeros(nf,nG);  % size: nf * nG

% function handle is not pde.Du
z0 = z1(1,:); v0 = fun(z0);
if length(v0)>1
    Du = fun;
    % nvec
    v12 = z2-z1; v13 = z3-z1;
    nvec = mycross(v12,v13,2);
    nvec = nvec./repmat(sqrt(sum(nvec.^2,2)),1,size(nvec,2));
    % Coef matrix: gN
    for p = 1:nG
        pz = lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
        Cmat(:,p) = sum(Du(pz).*nvec,2);
    end
    return;  % The remaining code will be neglected.
end

% function handle is not pde.Du
for p = 1:nG
    pz = lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
    Cmat(:,p) = fun(pz);
end

