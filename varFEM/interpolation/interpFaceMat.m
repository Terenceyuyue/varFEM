function Cmat = interpFaceMat(fun,Th,quadOrder)
%%interpFaceMat returns the matrix for numerical integration on boundary
% faces in 3-D
%
%  This function is designed to deal with terms containing normal
%  vector, for example, fun = pde.Du*n
%

if nargin==2, quadOrder = 3; end

node = Th.node;  
if isfield(Th,'on')
    on = Th.on;
    elem = Th.bdFaceType{on};
end
if isfield(Th,'elem'), elem = Th.elem; end  % provided by user
nf = size(elem,1); % 2-D mesh on th boundary

% Guass-Quadrature
[lambda,weight] = quadpts(quadOrder); ng = length(weight);

z1 = node(elem(:,1),:); 
z2 = node(elem(:,2),:); 
z3 = node(elem(:,3),:);
Cmat = zeros(nf,ng);  % size: nf * ng

% function handle is pde.Du
z0 = z1(1,:); v0 = fun(z0);
if length(v0)>1
    Du = fun;
    % nvec
    v12 = z2-z1; v13 = z3-z1;
    nvec = mycross(v12,v13,2);
    nvec = nvec./repmat(sqrt(sum(nvec.^2,2)),1,size(nvec,2));
    % Coef matrix: gN
    for p = 1:ng
        pz = lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
        Cmat(:,p) = sum(Du(pz).*nvec,2);
    end
    return;  % The remaining code will be neglected.
end

% function handle is not pde.Du
for p = 1:ng
    pz = lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
    Cmat(:,p) = fun(pz);
end