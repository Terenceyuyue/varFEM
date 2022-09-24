function [elemuhM,elemuhP] = elem2edgeInterp(wStr,Th,uh,Vh,quadOrder)

node = Th.node; elem = Th.elem;
NT = Th.NT; NE = Th.NE;

%% preparation for the computation
% basis functions
phi = Base2dBd(wStr,node,elem,Vh,quadOrder);
ng = size(phi{1},2)/3;
% auxiliary data
% totalEdge = sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2);
% [~, i1, totalJ] = unique(totalEdge,'rows');
% elem2edge = reshape(totalJ,NT,3);
% i2(totalJ) = 1:length(totalJ); i2 = i2(:); % second occurence by overlapping 
% bdEdgeIdx = (i1==i2);
elem2edge = Th.elem2edge;
bdEdgeIdx = Th.bdEdgeIdx;
% sgnelem
sgnelem = sign([elem(:,3)-elem(:,2), elem(:,1)-elem(:,3), elem(:,2)-elem(:,1)]);
E = false(NE,1); E(bdEdgeIdx) = 1; sgnbd = E(elem2edge);
sgnelem(sgnbd) = 0; % on the boundary of the domain
% elem2dof
elem2dof = dof2d(Th,Vh);

%% elementwise interior and exterior indices of quadrature points
elem2dofM = zeros(NT,3*ng); % M: minus -
%elem2dofP = zeros(NT,3*ng); % P: plus +
for i = 1:3  % i-th side
    ei = elem2edge(:,i); sgni = sgnelem(:,i);
    id = repmat(1:ng, NT,1);  
    id(sgni<0,:) = repmat((ng:-1:1)+ng*NE, sum(sgni<0), 1);
    elem2dofM(:,(1:ng)+(i-1)*ng) = id + (ei-1)*ng;
end
rep = ones(1,ng);
elem2dofsign = [sgnelem(:,rep),sgnelem(:,2*rep),sgnelem(:,3*rep)];
elem2dofP = elem2dofM + (-ng*NE)*(elem2dofsign<0) + ng*NE*(elem2dofsign>0); 

%% elementwise interior and exterior evaluations
elemuhM = zeros(NT,3*ng);
for p = 1:3*ng
    % interpolation at the p-th quadrture point
    for i = 1:length(phi)
        base = phi{i};
        elemuhM(:,p) = elemuhM(:,p) + uh(elem2dof(:,i)).*base(:,p);
    end
end
uhI = zeros(2*NE*ng,1);      
uhI(elem2dofM) = elemuhM;   
elemuhP = uhI(elem2dofP);  