function [elemuhP,elemuhM] = elem2edgeInterp(wStr,Th,uh,Vh,quadOrder)

node = Th.node; elem = Th.elem;
NT = Th.NT; NE = Th.NE;

%% preparation for the computation
% basis functions
phi = Base2dBd(wStr,node,elem,Vh,quadOrder);
ng = size(phi{1},2)/3;
% auxiliary data
totalEdge = sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2);
[~, i1, totalJ] = unique(totalEdge,'rows');
elem2edge = reshape(totalJ,NT,3);
i2(totalJ) = 1:length(totalJ); i2 = i2(:); % second occurence by overlapping 
bdEdgeIdx = (i1==i2);
% sgnelem
sgnelem = sign([elem(:,3)-elem(:,2), elem(:,1)-elem(:,3), elem(:,2)-elem(:,1)]);
E = false(NE,1); E(bdEdgeIdx) = 1; sgnbd = E(elem2edge);
sgnelem(sgnbd) = 1;
% elem2dof
elem2dof = dof2d(Th,Vh);

%% connectivity list of jump integral
elem2dofP = zeros(NT,3*ng); % P: plus +
%elem2dofM = zeros(NT,3*ng); % M: minus -
for i = 1:3  % i-th side
    ei = elem2edge(:,i); sgni = sgnelem(:,i);
    idP = repmat(1:ng, NT,1);  
    idP(sgni<0,:) = repmat((ng:-1:1)+NE*ng, sum(sgni<0), 1);
    elem2dofP(:,(1:ng)+(i-1)*ng) = idP + (ei-1)*ng;
end
index = (elem2dofP>ng*NE);
elem2dofM = elem2dofP + (-ng*NE)*index + ng*NE*(~index);

%% elementwise interpolant at quadrature points
elemuhP = zeros(NT,3*ng);
for p = 1:3*ng
    % interpolation at the p-th quadrture point
    for i = 1:length(phi)
        base = phi{i};
        elemuhP(:,p) = elemuhP(:,p) + uh(elem2dof(:,i)).*base(:,p);
    end
end
uhI = zeros(2*NE*ng,1);      
uhI(elem2dofP) = elemuhP;   
elemuhM = uhI(elem2dofM);  