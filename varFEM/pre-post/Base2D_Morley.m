% This is not a separate matlab script storing the base matrix of Morley element.
% see Base2D.m
%
% Copyright (C) Terence Yu.

%% structure
% numbers
N = size(node,1); NT = size(elem,1);
% edge,elem2edge,bdEdgeIdx
allEdge = uint32([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])]);
totalEdge = sort(allEdge,2);
[edge,i1,totalJ] = unique(totalEdge, 'rows');  NE = size(edge,1);
elem2edge  = reshape(totalJ, NT, 3);
% bdEdgeIdx
[~, i2] = unique(totalEdge(end:-1:1,:),'rows');
i2 = size(totalEdge,1)+1-i2;            % last occurrence
bdEdgeIdx = (i1==i2);

%% Preparation for the computation
% area
z1 = node(elem(:,1),:);
z2 = node(elem(:,2),:);
z3 = node(elem(:,3),:);
e1 = z2-z3; e2 = z3-z1; e3 = z1-z2; % ei = [xi, etai]
area = 0.5*(-e3(:,1).*e2(:,2)+e3(:,2).*e2(:,1));
% xi, eta
xi = [e1(:,1), e2(:,1), e3(:,1)];
eta = [e1(:,2), e2(:,2), e3(:,2)];
% sij
sb = zeros(NT,3,3);
for i = 1:3
    for j = 1:3
        sb(:,i,j) = -xi(:,i).*xi(:,j) - eta(:,i).*eta(:,j);
    end
end
% elementwise edge length
z1 = node(edge(:,1),:); z2 = node(edge(:,2),:);
he = sqrt(sum((z2-z1).^2,2)); L = he(elem2edge);
% elementwise sign of basis functions
sgnelem = sign([elem(:,3)-elem(:,2), elem(:,1)-elem(:,3), elem(:,2)-elem(:,1)]);
E = false(NE,1); E(bdEdgeIdx) = 1; sgnbd = E(elem2edge);
sgnelem(sgnbd) = 1;
sgnbase = ones(NT,6); sgnbase(:,4:6) = sgnelem;
% coefficients in the basis functions
ind = [1 2 3; 2 3 1; 3 1 2];  % rotation index
it = ind(:,1); jt = ind(:,2); kt = ind(:,3);
c0 = zeros(NT,3); c1 = c0; c2 = c0;
for i = 1:3
    j = ind(i,2); k = ind(i,3);
    c0(:,i) = 1./(2*area.^2);
    c1(:,i) = sb(:,i,j)./L(:,j).^2;
    c2(:,i) = sb(:,k,i)./L(:,k).^2;
end
c3 = c1+c2;

%% u.val
if mycontains(wStr,'.val')
    w1 = zeros(NT,nQuad); w2 = w1; w3 = w1; w4 = w1; w5 = w1; w6 = w1;
    w = {w1,w2,w3,w4,w5,w6};
    for p = 1:nQuad
        % basis functions at the p-th quadrture point
        base = zeros(NT,6);
        % i = 1,2,3
        base(:,it) = repmat(lambda(p,it).^2,NT,1) ...
            + c3.*repmat(lambda(p,jt).*lambda(p,kt),NT,1) ...
            + c2.*repmat(lambda(p,kt).*lambda(p,it),NT,1) ...
            + c1.*repmat(lambda(p,it).*lambda(p,jt),NT,1);
        % i = 4,5,6
        ci = 2*repmat(area,1,3)./L(:,it);
        base(:,3+it) = ci.*repmat(lambda(p,it).*(lambda(p,it)-1),NT,1);
        base = sgnbase.*base;  % equipped with global sign
        for i = 1:6
            w{i}(:,p) = base(:,i);
        end
    end
end

%% u.dxx
if mycontains(wStr,'.dxx')
    w1 = zeros(NT,nQuad); w2 = w1; w3 = w1; w4 = w1; w5 = w1; w6 = w1;
    w = {w1,w2,w3,w4,w5,w6};
    b11 = zeros(NT,6);
    % i = 1,2,3
    b11(:,it) = c0.* (eta(:,it).^2 - c1.*eta(:,jt).^2 - c2.*eta(:,kt).^2);
    % i = 4,5,6
    ci = 1./(repmat(area,1,3).*L(:,it));
    b11(:,3+it) = ci.*eta(:,it).^2;
    b11 = sgnbase.*b11;  % equipped with global sign
    for i = 1:6
        w{i} = repmat(b11(:,i),1,nQuad);
    end
end

%% u.dxy
if mycontains(wStr,'.dxy')
    w1 = zeros(NT,nQuad); w2 = w1; w3 = w1; w4 = w1; w5 = w1; w6 = w1;
    w = {w1,w2,w3,w4,w5,w6};
    b12 = zeros(NT,6);
    % i = 1,2,3
    b12(:,it) = -c0.* (xi(:,it).*eta(:,it) - c1.*xi(:,jt).*eta(:,jt) - c2.*xi(:,kt).*eta(:,kt));
    % i = 4,5,6
    ci = 1./(repmat(area,1,3).*L(:,it));
    b12(:,3+it) = -ci.*xi(:,it).*eta(:,it);
    b12 = sgnbase.*b12;  % equipped with global sign
    for i = 1:6
        w{i} = repmat(b12(:,i),1,nQuad);
    end
end

%% u.dyy
if mycontains(wStr,'.dyy')
    w1 = zeros(NT,nQuad); w2 = w1; w3 = w1; w4 = w1; w5 = w1; w6 = w1;
    w = {w1,w2,w3,w4,w5,w6};
    b22 = zeros(NT,6);
    % i = 1,2,3
    b22(:,it) = c0.* (xi(:,it).^2 - c1.*xi(:,jt).^2 - c2.*xi(:,kt).^2);
    % i = 4,5,6
    ci = 1./(repmat(area,1,3).*L(:,it));
    b22(:,3+it) = ci.*xi(:,it).^2;
    b22 = sgnbase.*b22;  % equipped with global sign
    for i = 1:6
        w{i} = repmat(b22(:,i),1,nQuad);
    end
end