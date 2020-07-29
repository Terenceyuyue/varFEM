function w = PlateBendingMorley1(node,elem,pde,bdStruct)

para = pde.para; f = pde.f; D = para.D; nu = para.nu;

%% elem2dof
% numbers
N = size(node,1); NT = size(elem,1); 
% edge and elem2edge
allEdge = uint32([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])]);
totalEdge = sort(allEdge,2);
[edge,~,totalJ] = unique(totalEdge, 'rows');
elem2edge  = reshape(totalJ, NT, 3);
% elem2dof
elem2dof = [elem, elem2edge+N];
% dof numbers
NE = size(edge,1);
Ndof = 6;  NNdof = N + NE;

% -------------- parameters used in the computation --------------------
auxG = auxgeometry(node,elem);  area = auxG.area;
x1 = node(elem(:,1),1); y1 = node(elem(:,1),2);
x2 = node(elem(:,2),1); y2 = node(elem(:,2),2);
x3 = node(elem(:,3),1); y3 = node(elem(:,3),2);
% xi, eta
xi = [x2-x3, x3-x1, x1-x2]; eta = [y2-y3, y3-y1, y1-y2];
% sij
sb = zeros(NT,3,3);
for i = 1:3
    j = 1:3;
    sb(:,i,j) = -xi(:,i).*xi(:,j) - eta(:,i).*eta(:,j);
end
% elementwise edge length
z1 = node(edge(:,1),:); z2 = node(edge(:,2),:);
he = sqrt(sum((z2-z1).^2,2)); L = he(elem2edge);
% signed elementwise edge length 
sgnelem = sign([elem(:,3)-elem(:,2), ...
                elem(:,1)-elem(:,3), ...
                elem(:,2)-elem(:,1)]);
bdEdgeIdx = bdStruct.bdEdgeIdx;
E = false(NE,1); E(bdEdgeIdx) = 1; sgnbd = E(elem2edge); 
sgnelem(sgnbd) = 1; 
sgnL = sgnelem.*L;
% coefficients in the basis functions
ind = [1 2 3; 2 3 1; 3 1 2];  % rotation index
it = ind(:,1); jt = ind(:,2); kt = ind(:,3);
c0 = zeros(NT,3); c1 = c0; c2 = c0;
for i = 1:3
    j = ind(i,2); k = ind(i,3);
    c0(:,i) = 1./(2*area.^2);
    c1(:,i) = sb(:,i,j)./sgnL(:,j).^2;
    c2(:,i) = sb(:,k,i)./sgnL(:,k).^2;
end
c3 = c1+c2;
% Gauss quadrature rule
[lambda,weight] = quadpts(2);
nQuad = length(weight);

%% Second derivatives of basis functions
b11 = zeros(NT,6); b22 = b11; b12 = b11; % [phi_i, i=1:6]
% i = 1,2,3
b11(:,it) = c0.* (eta(:,it).^2 - c1.*eta(:,jt).^2 - c2.*eta(:,kt).^2);
b22(:,it) = c0.* (xi(:,it).^2 - c1.*xi(:,jt).^2 - c2.*xi(:,kt).^2);
b12(:,it) = -c0.* (xi(:,it).*eta(:,it) - c1.*xi(:,jt).*eta(:,jt) - c2.*xi(:,kt).*eta(:,kt));
% i = 4,5,6
ci = 1./(area.*sgnL(:,it));
b11(:,3+it) = ci.*eta(:,it).^2;
b22(:,3+it) = ci.*xi(:,it).^2;
b12(:,3+it) = -ci.*xi(:,it).*eta(:,it);

%% First stiffness matrix
K = zeros(NT,Ndof^2);
for i = 1:Ndof
    j = 1:Ndof; jd = (i-1)*Ndof+1:i*Ndof;
    K(:,jd) = b11(:,i).*b11(:,j) + b22(:,i).*b22(:,j) ...
        + nu*(b11(:,i).*b22(:,j) + b22(:,i).*b11(:,j)) ...
        + 2*(1-nu)*b12(:,i).*b12(:,j);
end
K = D*area.*K;

%% Second stiffness matrix and load vector 
G = zeros(NT,Ndof^2); F = zeros(NT,Ndof);
if isnumeric(para.c), cf = @(xy) para.c+0*xy(:,1); end
for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    % basis functions at the p-th quadrture point
    base = zeros(NT,Ndof);
    % i = 1,2,3
    base(:,it) = lambda(p,it).^2 + c3.*lambda(p,jt).*lambda(p,kt) ...
        + c2.*lambda(p,kt).*lambda(p,it) + c1.*lambda(p,it).*lambda(p,jt);
    % i = 4,5,6
    ci = 2*area./sgnL(:,it);
    base(:,3+it) = ci.*lambda(p,it).*(lambda(p,it)-1);    
    % Second stiffness matrix
    for i = 1:Ndof
        j = 1:Ndof; jd = (i-1)*Ndof+1:i*Ndof;
        gs = cf(pxy).*base(:,i).*base(:,j);
        G(:,jd) = G(:,jd) + weight(p)*gs;
    end
    % load vector
    F = F + weight(p)*f(pxy).*base;
end
G = area.*G; F = area.*F;

%% Assemble stiffness matrix and load vector
ii = reshape(repmat(elem2dof, Ndof,1), [], 1);
jj = repmat(elem2dof(:), Ndof, 1);
kk = sparse(ii,jj,K(:)+G(:),NNdof,NNdof);
ff = accumarray(elem2dof(:), F(:), [NNdof 1]);

%% Dirichlet boundary conditions 
bdNodeIdx = bdStruct.bdNodeIdx; 
bdEdgeD = bdStruct.bdEdgeD;
g_D = pde.g_D;  Dw = pde.Dw; 

fixedNode = [bdNodeIdx; bdEdgeIdx+N];
isBdNode = false(N+NE,1); isBdNode(fixedNode) = true;
freeDof = ~isBdNode;

nodeD = node(bdNodeIdx,:); wD = g_D(nodeD); 
z1 = node(bdEdgeD(:,1),:); 
z2 = node(bdEdgeD(:,2),:); 
zc = (z1+z2)./2;
e = z1-z2;  % e = z2-z1
ne = [-e(:,2),e(:,1)];  ne = ne./he(bdEdgeIdx);  
wnD = sum(Dw(zc).*ne,2);
w = zeros(N+NE,1); w(fixedNode) = [wD;wnD];
ff = ff - kk*w;

%% Solver
w(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
