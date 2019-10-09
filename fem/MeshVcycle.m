function [node,elem,Pro,Res] = MeshVcycle(node0,elem0,J)

%J = 4; % level length, J>=2
Pro = cell(J-1,1); Res = cell(J-1,1);
for j = 2:J
    % node and elem in each level
    node1 = node0; elem1 = elem0;
    xc = (node1(elem1(:,1))+node1(elem1(:,2)))/2;
    N = size(node1,1); nel = size(elem1,1);  idc = (1:nel)' + N;
    node = [node1;xc];
    elem = [[elem1(:,1),idc];    % left subintervals
           [idc, elem1(:,2)]];   % right subintervals
    node0 = node; elem0 = elem;
    
    % HB in the current level
    HB = zeros(nel,3);
    HB(:,1) = idc; HB(:,2:3) = elem1;
    
    % Prolongation matrix
    Coarse2Fine = (1:N)'; CoarseId = (1:N)';
    nCoarse = N; nf = nel; nFine = nCoarse+nf;    
    ii = [Coarse2Fine; HB(:,1); HB(:,1)];
    jj = [CoarseId; HB(:,2); HB(:,3)];
    ss = [ones(nCoarse,1); 0.5*ones(nf,1); 0.5*ones(nf,1)];
    Pro{j-1} = sparse(ii,jj,ss,nFine,nCoarse);
    Res{j-1} = Pro{j-1}';
end
