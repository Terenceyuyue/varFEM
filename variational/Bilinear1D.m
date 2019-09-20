function kk = Bilinear1D(Th,varargin)

if nargin==3
    Trial = varargin{1}; Test = varargin{2};
    Coef = num2cell(1:length(Trial));
elseif nargin==4
    Coef = varargin{1}; Trial = varargin{2}; Test = varargin{3};
else
    error('Check the input');
end

node = Th.node; elem = Th.elem;
N = size(node,1); nel = size(elem,1); Ndof = 2;
xa = node(elem(:,1)); xb = node(elem(:,2)); 

K = zeros(nel,2*Ndof); % straighten for easy assembly
J = Jacobian1D(xa,xb);  % Jacobian on all elements
[r, w] = GaussQuad(4,0,1); % on [0,1];
xx = xa+(xb-xa)*r';  % quadrature points on all elements
for s = 1:length(Coef)  % s-th component of bilinear form
    u = Trial{s}; v = Test{s};
    % Trial and Test functions
    [u1,u2] = Refbase1D(u); [v1,v2] = Refbase1D(v);
    % Stiffness matrix on reference interval [0,1]
    cf = Coef{s};  
    if isa(cf,'double'), cf = @(x) cf+0*x; end
    cf = cf(xx);
    v1 = v1(r); v2 = v2(r); u1 = u1(r); u2 = u2(r);
    ke = w.*[v1.*u1, v1.*u2, v2.*u1, v2.*u2];
    Kref = cf*ke;
    % Stiffness matrix on all elements    
    cu = Trans1D(u,xa,xb);  cv = Trans1D(v,xa,xb); 
    K = K + cv.*cu.*J.*Kref; 
end

% sparse indices
nnz = nel*Ndof^2;
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = K(:);
id = 0;
for i = 1:Ndof
    for j = 1:Ndof
        ii(id+1:id+nel) = elem(:,i);   % zi
        jj(id+1:id+nel) = elem(:,j);   % zj
        id = id + nel;
    end
end

% stiffness matrix
kk = sparse(ii,jj,ss,N,N);
