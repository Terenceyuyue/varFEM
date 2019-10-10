clc;clear;close all

% -------- Mesh and boudary conditions --------
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
Nx = 10; Ny = 10; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);
bdStruct = setboundary(node,elem);

% ----------------- PDE data -------------------
pde = PlateBendingData;
para = pde.para;
f = pde.f;

D = para.D; nu = para.nu;
R = D*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];

% ------------- Sparse assembling indices ------------
N = size(node,1); NT = size(elem,1); Ndof = 9;
elem2 = zeros(NT,Ndof);
elem2(:,[1,4,7]) = elem;
elem2(:,[2,5,8]) = elem + N;
elem2(:,[3,6,9]) = elem + 2*N;
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
% xi, eta
x1 = node(elem(:,1),1); y1 = node(elem(:,1),2);
x2 = node(elem(:,2),1); y2 = node(elem(:,2),2);
x3 = node(elem(:,3),1); y3 = node(elem(:,3),2);
xi = [x2-x3, x3-x1, x1-x2]; eta = [y2-y3, y3-y1, y1-y2];
% Gauss quadrature rule
n = 3;
[lambda,weight] = quadpts(n);
nI = length(weight); % number of integration points

% ------- Compute second derivatives of triple homogeneous polynomials ----------
% Derivative of \lambda_i
[Dphi,area] = gradbasis(node,elem);


K = zeros(NT,Ndof^2);
for iel = 1:NT
    Dqua = @(i,j) Dquad(i,j,iel,Dphi,lambda);
    Dtri = @(i,j,k) Dtriple(i,j,k,iel,Dphi,lambda);    
    
    Dbase1 = cell(3,1); Dbase2 = cell(3,1); Dbase3 = cell(3,1);
    for i = 1:3
        % \phi (11,22,12)
        Dbase1{i} = 3*Dqua(i,i)-2*Dtri(i,i,i)+2*Dtri(1,2,3);
        
        % \psi (11,22,12)
        c1 = xi(iel,1)*ones(1,3);
        c2 = xi(iel,2)*ones(1,3);
        c3 = [xi(iel,2),xi(iel,1),0];
        c4 = [0.5*xi(iel,1)+xi(iel,2),  -(0.5*xi(iel,2)+xi(iel,1)),  0.5*(xi(iel,1)-xi(iel,2))];
        Dbase2{i} = c1(i)*Dtri(i,i,2)-c2(i)*Dtri(i,i,1)-c3(i)*Dqua(i,i)...
            + c4(i)*Dtri(1,2,3);
        
        % \zeta (11,22,12)
        c1 = eta(iel,1)*ones(1,3);
        c2 = eta(iel,2)*ones(1,3);
        c3 = [eta(iel,2),eta(iel,1),0];
        c4 = [0.5*eta(iel,1)+eta(iel,2), -(0.5*eta(iel,2)+eta(iel,1)), 0.5*(eta(iel,1)-eta(iel,2))];        
        Dbase3{i} = c1(i)*Dtri(i,i,2)-c2(i)*Dtri(i,i,1)-c3(i)*Dqua(i,i)...
            + c4(i)*Dtri(1,2,3);
    end
    
    Ke = cell(3,3);
    for s = 1:3
        for t = 1:3
            Bs = zeros(3,3); Bt = zeros(3,3); Kst = zeros(3,3);
            for m = 1:nI  % »ý·ÖµãÐòºÅ
                Bs(:,1) = -Dbase1{s}(m,:);  Bt(:,1) = -Dbase1{t}(m,:);
                Bs(:,2) = -Dbase2{s}(m,:);  Bt(:,2) = -Dbase2{t}(m,:);
                Bs(:,3) = -Dbase3{s}(m,:);  Bt(:,3) = -Dbase3{t}(m,:);  
                Bs(3,:) = 2*Bs(3,:);   Bt(3,:) = 2*Bt(3,:);
                Kst = Kst + weight(m)*Bs'*R*Bt;
            end
            Ke{s,t} = area(iel)*Kst;
        end
    end
    Ke = cell2mat(Ke);
    K(iel,:) = reshape(Ke',1,[]);
end
kk = sparse(ii,jj,K(:),3*N,3*N);


% -------------- ÓÒ¶Ë -------------
F = zeros(NT,Ndof);
base = zeros(nI,Ndof);
for iel = 1:NT
    vK = node(elem(iel,:),:); % vertices of K
    xy = lambda*vK;
    fxy = f(xy);
    
    for i = 1:3
        % \phi
        base(:,i) = lambda(:,i).^2.*(3-2*lambda(:,i)) + 2*lambda(:,1).*lambda(:,2).*lambda(:,3);
        % \psi
        c1 = xi(iel,1)*ones(1,3);
        c2 = xi(iel,2)*ones(1,3);
        c3 = [xi(iel,2),xi(iel,1),0];
        c4 = [0.5*xi(iel,1)+xi(iel,2),-(0.5*xi(iel,2)+xi(iel,1)),0.5*(xi(iel,1)-xi(iel,2))];
        base(:,3+i) = lambda(:,i).^2.*(c1(i)*lambda(:,2)-c2(i)*lambda(:,1)-c3(i)) ...
            + c4(i)*lambda(:,1).*lambda(:,2).*lambda(:,3);
        % \zeta
        c1 = eta(iel,1)*ones(1,3);
        c2 = eta(iel,2)*ones(1,3);
        c3 = [eta(iel,2),eta(iel,1),0];
        c4 = [0.5*eta(iel,1)+eta(iel,2),-(0.5*eta(iel,2)+eta(iel,1)),0.5*(eta(iel,1)-eta(iel,2))];
        base(:,6+i) = lambda(:,i).^2.*(c1(i)*lambda(:,2)-c2(i)*lambda(:,1)-c3(i)) ...
            + c4(i)*lambda(:,1).*lambda(:,2).*lambda(:,3);
    end    
    base = base(:,[1 4 7 2 5 8 3 6 9]);
    fv = fxy.*base;
    F(iel,:) = area(iel)*weight*fv;
end
ff = accumarray(elem2(:), F(:), [3*N 1]);

% ------------ Dirichlet boundary conditions ----------------
eD = bdStruct.eD;
id = [eD; eD+N; eD+2*N];
isBdDof = false(3*N,1); isBdDof(id) = true;
bdDof = find(isBdDof); freeDof = find(~isBdDof);

g_D = pde.g_D;  pD = node(eD,:);  Dw = pde.Dw;
wD = g_D(pD);  Dw = Dw(pD); wxD = Dw(:,1); wyD = Dw(:,2);

w = zeros(3*N,1); w(bdDof) = [wD;wxD;wyD];

ff = ff - kk*w;

% ------------------ Solver -------------------
w(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
w = w(1:N);


% --------- error analysis -------
uexact = pde.uexact;  we = uexact(node);
figure,
subplot(1,2,1), showsolution(node,elem,w);
zlabel('w');
subplot(1,2,2), showsolution(node,elem,we);
zlabel('we');
Eabs = w-we;  % Absolute errors
figure,showsolution(node,elem,Eabs); zlim('auto');

format shorte
Err = norm(Eabs)./norm(we)

