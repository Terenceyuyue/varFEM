function w = PlateBendingAdini(node,elem,pde,bdStruct)

para = pde.para; f = pde.f; D = para.D; nu = para.nu;
% -------- Sparse assembling indices -----------
N = size(node,1); NT = size(elem,1); Ndof = 4*3;
elem2 = [elem, elem+N, elem + 2*N];

nnz = NT*Ndof^2;
ii = zeros(nnz,1); jj = zeros(nnz,1);
id = 0;
for i = 1:Ndof
    for j = 1:Ndof
        ii(id+1:id+NT) = elem2(:,i);
        jj(id+1:id+NT) = elem2(:,j);
        id = id + NT;
    end
end

% ---------------- parameters used in computation -----------------------
xiv = [-1 1 1 -1]; etav = [-1 -1 1 1]; % vertices of [-1,1]^2
x1 = node(elem(:,1),1); y1 = node(elem(:,1),2);
x2 = node(elem(:,2),1); y4 = node(elem(:,4),2);
a = x2-x1;  b = y4-y1;
% quadature points
[xx,ww] = quadpts1(5);  nr = length(ww); nI = nr^2;
xi = reshape(repmat(xx',nr,1),1,[]); 
eta = reshape(ones(nr,1)*xx,1,[]);
weight = reshape(ww'*ww,1,[]); 


% ------- second derivatives of basis functions ----------
b11 = zeros(NT,Ndof,nI); b22 = b11; b12 = b11;
for i = 1:4
    % \phi
    b11(:,i,:) = -3./a.^2*xiv(i)*xi.*(1+etav(i)*eta);
    b22(:,i,:) = -3./b.^2*etav(i)*eta.*(1+xiv(i)*xi);
    b12(:,i,:) = 1./(2*a.*b)*xiv(i)*etav(i)*(4-3*(xi.^2+eta.^2));
    % \psi
    b11(:,4+i,:) = 1./(2*a)*(1+etav(i)*eta).*(xiv(i)+3*xi);
    b22(:,4+i,:) = 0;
    b12(:,4+i,:) = -1./(4*b)*etav(i)*(1-2*xiv(i)*xi-3*xi.^2);
    % \zeta
    b11(:,8+i,:) = 0;
    b22(:,8+i,:) = 1./(2*b)*(1+xiv(i)*xi).*(etav(i)+3*eta);
    b12(:,8+i,:) = -1./(4*a)*xiv(i)*(1-2*etav(i)*eta-3*eta.^2);
end

% ----------- First stiffness matrix -----------
K = zeros(NT,Ndof^2);
for i = 1:Ndof
    j = 1:Ndof; jd = (i-1)*Ndof+1:i*Ndof;
    for p = 1:nI
        Kp = b11(:,i,p).*b11(:,j,p) + b22(:,i,p).*b22(:,j,p) ...
            + nu*(b11(:,i,p).*b22(:,j,p) + b22(:,i,p).*b11(:,j,p)) ...
            + 2*(1-nu)*b12(:,i,p).*b12(:,j,p);
        K(:,jd) = K(:,jd) + weight(p)*Kp;
    end
end
K = D*(a.*b)./4.*K;

% ------------- Second stiffness matrix and load vector ------------
G = zeros(NT,Ndof^2); F = zeros(NT,Ndof); 
if isnumeric(para.c), cf = @(xy) para.c+0*xy(:,1); end
x0 = (x1+x2)./2; y0 = (y1+y4)./2;
for p = 1:nI
    % quadrature points in the x-y coordinate
    x = a*xi(p)/2 + x0; y = b*eta(p)/2 + y0;
    pxy = [x,y];
    % basis functions at the p-th quadrture point
    base = zeros(NT,Ndof); 
    for i = 1:4
        tp = (1+xiv(i)*xi(p))*(1+etav(i)*eta(p));
        % \phi
        base(:,i) = 1/8*tp*(2+xiv(i)*xi(p)+etav(i)*eta(p)-xi(p)^2-eta(p)^2);
        % \psi
        base(:,4+i) = -1/16*a*xiv(i)*tp*(1-xi(p)^2);
        % \zeta
        base(:,8+i) = -1/16*b*etav(i)*tp*(1-eta(p)^2);
     end    
    % Second stiffness matrix
    for i = 1:Ndof
        j = 1:Ndof; jd = (i-1)*Ndof+1:i*Ndof;
        gs = cf(pxy).*base(:,i).*base(:,j);
        G(:,jd) = G(:,jd) + weight(p)*gs;
    end    
    % load vector
    F = F + weight(p)*f(pxy).*base;    
end
G = (a.*b)./4.*G; F = (a.*b)./4.*F;

% ------------- Assemble stiffness matrix and load vector ------------
kk = sparse(ii,jj,K(:)+G(:),3*N,3*N);
ff = accumarray(elem2(:), F(:), [3*N 1]);

% ------------ Dirichlet boundary conditions ----------------
eD = bdStruct.eD;
id = [eD; eD+N; eD+2*N];
isBdDof = false(3*N,1); isBdDof(id) = true;
bdDof = isBdDof; freeDof = (~isBdDof);
g_D = pde.g_D;  pD = node(eD,:);  Dw = pde.Dw;
wD = g_D(pD);  Dw = Dw(pD); %wxD = Dw(:,1); wyD = Dw(:,2);
w = zeros(3*N,1); w(bdDof) = [wD;Dw(:)];
ff = ff - kk*w;

% ------------------ Solver -------------------
w(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
w = w(1:N);