function [node,elem,Pro,Res] = Mesh1DVcycle(node,elem,J)

if J<=1
    Pro = []; Res = []; return;
end

Pro = cell(J-1,1); Res = cell(J-1,1);
for j = 2:J
    
    N = size(node,1); NT = size(elem,1); Ndof = 2;
    % HB in the current level    
    HB = zeros(NT,3);
    HB(:,1) = (1:NT)'+N;  HB(:,2:3) = elem;
    
    % Prolongation matrix
    Coarse2Fine = (1:N)'; CoarseId = (1:N)';
    nCoarse = N; nf = NT; nFine = nCoarse+nf;
    ii = [Coarse2Fine; HB(:,1); HB(:,1)];
    jj = [CoarseId; HB(:,2); HB(:,3)];
    ss = [ones(nCoarse,1); 0.5*ones(nf,1); 0.5*ones(nf,1)];
    Pro{j-1} = sparse(ii,jj,ss,nFine,nCoarse);
    Res{j-1} = Pro{j-1}';
    
    % add new nodes
    node(N+1:N+NT) = (node(elem(:,1))+node(elem(:,2)))/2;
    
    % add new elements
    % 1 --- 3 --- 2
    t = 1:NT; p = zeros(NT,2*Ndof);
    p(:,1:2) = elem;  p(:,3) = (1:NT)' + N;
    elem(t,:) = [p(t,1), p(t,3)];
    elem(NT+1:2*NT,:) = [p(t,3), p(t,2)];     
end

