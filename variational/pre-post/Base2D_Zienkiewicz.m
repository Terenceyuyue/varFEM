% This is not a separate matlab script storing the base matrix of Zienkiewicz element.
% see Base2D.m
%
% Copyright (C) Terence Yu.

%% parameters used in computation
% xi, eta
x1 = node(elem(:,1),1); y1 = node(elem(:,1),2);
x2 = node(elem(:,2),1); y2 = node(elem(:,2),2);
x3 = node(elem(:,3),1); y3 = node(elem(:,3),2);
xi = [x2-x3, x3-x1, x1-x2]; eta = [y2-y3, y3-y1, y1-y2];
% coefficients in the basis functions
c1 = xi(:,1);   d1 = eta(:,1); 
c2 = -xi(:,2);  d2 = -eta(:,2); 
c3 = -[c2,c1,zeros(NT,1)]; d3 = -[d2,d1,zeros(NT,1)];
c4 = [0.5*c1-c2,  0.5*c2-c1,  0.5*(c1+c2)];
d4 = [0.5*d1-d2,  0.5*d2-d1,  0.5*(d1+d2)];

%% second derivatives of homogeneous polynomials
Dquad = @(i,j) quadbasis(i,j,Dlambda,lambda);
Dtri = @(i,j,k) tribasis(i,j,k,Dlambda,lambda);

%% u.val
if mycontains(wStr,'.val')
    w1 = zeros(NT,nQuad); w2 = w1; w3 = w1;
    w4 = w1; w5 = w1; w6 = w1; w7 = w1; w8 = w1; w9 = w1;
    w = {w1,w2,w3,w4,w5,w7,w8,w9};
    for p = 1:nQuad
        % basis functions at the p-th quadrture point
        % basis functions at the p-th quadrture point
        base = zeros(NT,9);
        for i = 1:3
            % \phi
            lam123 = lambda(p,1).*lambda(p,2).*lambda(p,3);
            base(:,i) = lambda(p,i).^2.*(3-2*lambda(p,i)) + 2*lam123;
            % \psi
            base(:,3+i) = lambda(p,i).^2.*(c1*lambda(p,2) + c2*lambda(p,1) + c3(:,i)) ...
                + c4(:,i)*lam123;
            % \zeta
            base(:,6+i) = lambda(p,i).^2.*(d1*lambda(p,2) + d2*lambda(p,1) + d3(:,i)) ...
                + d4(:,i)*lam123;
        end
        for i = 1:9
            w{i}(:,p) = base(:,i);
        end
    end
end

%% u.dxx
if mycontains(wStr,'.dxx')
    w1 = zeros(NT,nQuad); w2 = w1; w3 = w1; 
    w4 = w1; w5 = w1; w6 = w1; w7 = w1; w8 = w1; w9 = w1;
    w = {w1,w2,w3,w4,w5,w7,w8,w9};    
    for i = 1:3
        % \phi 
        DD = 3*Dquad(i,i)-2*Dtri(i,i,i)+2*Dtri(1,2,3);
        w{i} = DD(:,:,1);
        % \psi
        cc1(:,:,1) = repmat(c1,1,nQuad);      cc1(:,:,2) = cc1(:,:,1); cc1(:,:,3) = cc1(:,:,1);
        cc2(:,:,1) = repmat(c2,1,nQuad);      cc2(:,:,2) = cc2(:,:,1); cc2(:,:,3) = cc2(:,:,1);
        cc3(:,:,1) = repmat(c3(:,i),1,nQuad); cc3(:,:,2) = cc3(:,:,1); cc3(:,:,3) = cc3(:,:,1);
        cc4(:,:,1) = repmat(c4(:,i),1,nQuad); cc4(:,:,2) = cc4(:,:,1); cc4(:,:,3) = cc4(:,:,1);
        DD = cc1.*Dtri(i,i,2)+ cc2.*Dtri(i,i,1) + cc3.*Dquad(i,i) + cc4.*Dtri(1,2,3);
        w{3+i} = DD(:,:,1);
        % \zeta
        dd1(:,:,1) = repmat(d1,1,nQuad);      dd1(:,:,2) = dd1(:,:,1); dd1(:,:,3) = dd1(:,:,1);
        dd2(:,:,1) = repmat(d2,1,nQuad);      dd2(:,:,2) = dd2(:,:,1); dd2(:,:,3) = dd2(:,:,1);
        dd3(:,:,1) = repmat(d3(:,i),1,nQuad); dd3(:,:,2) = dd3(:,:,1); dd3(:,:,3) = dd3(:,:,1);
        dd4(:,:,1) = repmat(d4(:,i),1,nQuad); dd4(:,:,2) = dd4(:,:,1); dd4(:,:,3) = dd4(:,:,1);
        DD = dd1.*Dtri(i,i,2) + dd2.*Dtri(i,i,1) + dd3.*Dquad(i,i) + dd4.*Dtri(1,2,3);
        w{6+i} = DD(:,:,1);
    end
end

%% u.dyy
if mycontains(wStr,'.dyy')
    w1 = zeros(NT,nQuad); w2 = w1; w3 = w1; 
    w4 = w1; w5 = w1; w6 = w1; w7 = w1; w8 = w1; w9 = w1;
    w = {w1,w2,w3,w4,w5,w7,w8,w9};    
    for i = 1:3
        % \phi 
        DD = 3*Dquad(i,i)-2*Dtri(i,i,i)+2*Dtri(1,2,3);
        w{i} = DD(:,:,2);
        % \psi
        cc1(:,:,1) = repmat(c1,1,nQuad);      cc1(:,:,2) = cc1(:,:,1); cc1(:,:,3) = cc1(:,:,1);
        cc2(:,:,1) = repmat(c2,1,nQuad);      cc2(:,:,2) = cc2(:,:,1); cc2(:,:,3) = cc2(:,:,1);
        cc3(:,:,1) = repmat(c3(:,i),1,nQuad); cc3(:,:,2) = cc3(:,:,1); cc3(:,:,3) = cc3(:,:,1);
        cc4(:,:,1) = repmat(c4(:,i),1,nQuad); cc4(:,:,2) = cc4(:,:,1); cc4(:,:,3) = cc4(:,:,1);
        DD = cc1.*Dtri(i,i,2)+ cc2.*Dtri(i,i,1) + cc3.*Dquad(i,i) + cc4.*Dtri(1,2,3);
        w{3+i} = DD(:,:,2);
        % \zeta
        dd1(:,:,1) = repmat(d1,1,nQuad);      dd1(:,:,2) = dd1(:,:,1); dd1(:,:,3) = dd1(:,:,1);
        dd2(:,:,1) = repmat(d2,1,nQuad);      dd2(:,:,2) = dd2(:,:,1); dd2(:,:,3) = dd2(:,:,1);
        dd3(:,:,1) = repmat(d3(:,i),1,nQuad); dd3(:,:,2) = dd3(:,:,1); dd3(:,:,3) = dd3(:,:,1);
        dd4(:,:,1) = repmat(d4(:,i),1,nQuad); dd4(:,:,2) = dd4(:,:,1); dd4(:,:,3) = dd4(:,:,1);
        DD = dd1.*Dtri(i,i,2) + dd2.*Dtri(i,i,1) + dd3.*Dquad(i,i) + dd4.*Dtri(1,2,3);
        w{6+i} = DD(:,:,2);
    end
end

%% u.dxy
if mycontains(wStr,'.dxy')
    w1 = zeros(NT,nQuad); w2 = w1; w3 = w1; 
    w4 = w1; w5 = w1; w6 = w1; w7 = w1; w8 = w1; w9 = w1;
    w = {w1,w2,w3,w4,w5,w7,w8,w9};    
    for i = 1:3
        % \phi 
        DD = 3*Dquad(i,i)-2*Dtri(i,i,i)+2*Dtri(1,2,3);
        w{i} = DD(:,:,3);
        % \psi
        cc1(:,:,1) = repmat(c1,1,nQuad);      cc1(:,:,2) = cc1(:,:,1); cc1(:,:,3) = cc1(:,:,1);
        cc2(:,:,1) = repmat(c2,1,nQuad);      cc2(:,:,2) = cc2(:,:,1); cc2(:,:,3) = cc2(:,:,1);
        cc3(:,:,1) = repmat(c3(:,i),1,nQuad); cc3(:,:,2) = cc3(:,:,1); cc3(:,:,3) = cc3(:,:,1);
        cc4(:,:,1) = repmat(c4(:,i),1,nQuad); cc4(:,:,2) = cc4(:,:,1); cc4(:,:,3) = cc4(:,:,1);
        DD = cc1.*Dtri(i,i,2)+ cc2.*Dtri(i,i,1) + cc3.*Dquad(i,i) + cc4.*Dtri(1,2,3);
        w{3+i} = DD(:,:,3);
        % \zeta
        dd1(:,:,1) = repmat(d1,1,nQuad);      dd1(:,:,2) = dd1(:,:,1); dd1(:,:,3) = dd1(:,:,1);
        dd2(:,:,1) = repmat(d2,1,nQuad);      dd2(:,:,2) = dd2(:,:,1); dd2(:,:,3) = dd2(:,:,1);
        dd3(:,:,1) = repmat(d3(:,i),1,nQuad); dd3(:,:,2) = dd3(:,:,1); dd3(:,:,3) = dd3(:,:,1);
        dd4(:,:,1) = repmat(d4(:,i),1,nQuad); dd4(:,:,2) = dd4(:,:,1); dd4(:,:,3) = dd4(:,:,1);
        DD = dd1.*Dtri(i,i,2) + dd2.*Dtri(i,i,1) + dd3.*Dquad(i,i) + dd4.*Dtri(1,2,3);
        w{6+i} = DD(:,:,3);
    end
end