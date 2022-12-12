function [elemuhM,elemuhP,elemQuadnx,elemQuadny] = elem2edgeInterp(wStr,Th,uh,Vh,quadOrder)

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
sgnelem = sign(elem(:,[3,1,2])-elem(:,[2,3,1]));
E = false(NE,1); E(bdEdgeIdx) = 1; sgnbd = E(elem2edge);
sgnelem(sgnbd) = 0; % on the boundary of the domain
% elem2dof
elem2dof = dof2d(Th,Vh);

%% elementwise interior and exterior indices of quadrature points
elemQuadM = zeros(NT,3*ng); % M: minus -
%elem2dofP = zeros(NT,3*ng); % P: plus +
for i = 1:3  % i-th side
    ei = elem2edge(:,i); sgni = sgnelem(:,i);
    id = repmat(1:ng, NT,1);  
    id(sgni<0,:) = repmat((ng:-1:1)+ng*NE, sum(sgni<0), 1);
    elemQuadM(:,(1:ng)+(i-1)*ng) = id + (ei-1)*ng;
end
elemQuadP = elemQuadM + (-ng*NE)*(elemQuadM>ng*NE) + ng*NE*(elemQuadM<=ng*NE); 

%% elementwise interior and exterior evaluations
% interior evaluations
elemuhM = zeros(NT,3*ng);
for p = 1:3*ng
    % interpolation at the p-th quadrture point
    for i = 1:length(phi)
        if sum(uh(elem2dof(:,i)),"all") == 0, continue; end  
        base = phi{i};
        elemuhM(:,p) = elemuhM(:,p) + uh(elem2dof(:,i)).*base(:,p);
    end
end
% exterior evaluations
uhI = zeros(2*NE*ng,1);      
uhI(elemQuadM) = elemuhM;  % the exterior values on the boundary edges remain zero  
elemuhP = uhI(elemQuadP);  

%% elementwise unit normal vectors of edges on quadrature points
if nargout==4
    rep = ones(1,ng);
    z1 = node(elem(:,1),:);  
    z2 = node(elem(:,2),:); 
    z3 = node(elem(:,3),:);
    e1 = z2-z3;   e1 = e1./vecnorm(e1,2,2); % -e1
    e2 = z3-z1;   e2 = e2./vecnorm(e2,2,2); 
    e3 = z1-z2;   e3 = e3./vecnorm(e3,2,2); 
    elemQuadnx = -[e1(:,2*rep),e2(:,2*rep),e3(:,2*rep)];
    elemQuadny = [e1(:,rep),e2(:,rep),e3(:,rep)];
end