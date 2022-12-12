function [edgeBaseM,edgeBaseP,edge2dof] = edgeBase(wStr,Th,Vh,quadOrder)
% This function will greatly simplify the computation for bilinear forms
% involving jump and average terms.
% We compute and assemble the stiffness matrix as for bilinear forms defined on elements.
% This is based on the following correspondences:
%  -  edge e  --->  element K
%  -  basis functions on K1 \cup K2   --->   basis functions on K,
% where K1 and K2 share e as an edge.
%
% Therefore, we further introduce the following correspondences:
%  - connectivity: edge2dof  --->  elem2dof
%  - Basis matrix for numerical integration: edgeBase --->  Base2d
% Here, edgeBase includes the left and right evaluations, indicating by M
% (minus) and P (Plus), respectively.
%
% Then the bilinear forms can be assembled as in assem2d.m, see main_biharmonic_C0IP.m
% for example, where the the symmetric quadratic C0 interior penalty method for the
% biharmonic equation: \Delta^2 u = f
%

node = Th.node; elem = Th.elem;
if ~iscell(wStr), wStr = {wStr}; end
nwStr = length(wStr);

%% Preparation for the computation
% numbers
Ndofe = 9; % 9 dofs for the two elements sharing e as an edge
NE = Th.NE;
[~,w1d] = quadpts1(quadOrder);  ng = length(w1d);
% aux
edge2elem = Th.edge2elem;
bdEdgeIdx = Th.bdEdgeIdx;
% sgnelem
sgnelem = sign(elem(:,[3 1 2])-elem(:,[2 3 1]));
E = false(NE,1); E(bdEdgeIdx) = 1; sgnbd = E(Th.elem2edge);
sgnelem(sgnbd) = 0; % on the domain boundary
% elem2dof
[elem2dof,Ndof] = dof2d(Th,Vh);
% false for the boundary edges
isNotBdEdge = true(NE,1);  isNotBdEdge(bdEdgeIdx) = false;

%% Information of left and right elements
% left and right elements
k1 = edge2elem(:,1);  k2 = edge2elem(:,2);
% local index of the edge on left or right elements
e1 = edge2elem(:,3);   e2 = edge2elem(:,4);
% global dofs on left and right elements
indexk1 = elem2dof(k1,:);   indexk2 = elem2dof(k2,:);
% sign of the edge on the left element
sgnk1 = subMat(sgnelem,k1,e1); % subMat(A,i,j) gets the (i,j)-th entry

%% Macro indices and Data structure edge2dof
edge2dof = zeros(NE,Ndofe);  % similar to elem2dof
idMacro = zeros(NE, 2*Ndof); % macro indices of left and right basis functions
for ie = 1:NE
    [index,~, Ia] = unique([indexk1(ie,:),indexk2(ie,:)]);
    idMacro(ie, :) = Ia;
    index = index(:)';  Nv = length(index);
    index1 = index;
    if Nv<Ndofe
        index1 = [index, ones(1,Ndofe-Nv)]; % 1：virtual dofs for boundary edges
    end
    edge2dof(ie, :) = index1;    
end

%% Row and column subscripts used in subMat
% left: base(k1, (1:ng)+(e1-1)*ng)
Lrows = repmat(k1,1,ng);
Lcols = repmat(1:ng, NE, 1) + repmat((e1-1)*ng,1,ng);
% right: base(k2, (1:ng)+(e2-1)*ng);
Rrows = repmat(k2,1,ng);
Rcols = repmat(1:ng, NE, 1) + repmat((e2-1)*ng,1,ng);
% reverse order 
ReverseRows = repmat((1:Ndofe*NE)', 1, 2*ng);
ReverseCols = repmat(sgnk1>=0,1,2*ng).*[(1:ng), (ng:-1:1)+ng] ...
    + repmat(sgnk1<0,1,2*ng).*[(1:ng)+ng, (ng:-1:1)];  % (NE, 2*ng)
ReverseCols = kron(ReverseCols, ones(Ndofe,1));
% increments for left and right dofs
rP = (0: Ndofe: (NE-1)*Ndofe)';

%% Basis matrix for numerical integration on edges
% edgePhi has NE row blocks and the ib-th block is for the ib-th basis
% function. Let Phi be some row block. Then it is of size (9, 2*ng), where
% 9 is the number of the dofs on the macro element . The first ng columns
% are for the left element and the remainings are for the right element.
edgePhi = zeros(Ndofe*NE, 2*ng, nwStr); 
for ss = 1:nwStr 

    % basis functions for the ss-th string
    phi = Base2dBd(wStr{ss},node,elem,Vh,quadOrder); 

    % loop of basis functions on the left or right elements
    for ib = 1:Ndof 
        % the current basis functions
        base = phi{ib}; % (NT, 3*ng)
        % indices in edgePhi     
        rowL = idMacro(:, ib) + rP;      
        rowR = idMacro(:, ib+Ndof) + rP; 
        % left: base(k1, (1:ng)+(e1-1)*ng)        
        edgePhi(rowL, 1:ng, ss) = subMat(base, Lrows, Lcols); 
        % right: base(k2, (1:ng)+(e2-1)*ng)        
        edgePhi(rowR(isNotBdEdge), (1:ng)+ng, ss) = ...  % k1~=k2
            subMat(base, Rrows(isNotBdEdge,:), Rcols(isNotBdEdge,:));
    end

    % reorder w.r.t the positive sign
    edgePhi(:,:,ss) = subMat(edgePhi(:,:,ss), ReverseRows, ReverseCols);
end

%% Stored in cell form
edgeBaseM = cell(1,nwStr);  edgeBaseP = cell(1,nwStr); 
for ss = 1:nwStr
    phiM = cell(1,Ndofe);  phiP = cell(1,Ndofe);
    for i = 1:Ndofe
        phiM{i} = edgePhi(i:Ndofe:end, 1:ng, ss);
        phiP{i} = edgePhi(i:Ndofe:end, (1:ng)+ng, ss);
    end
    edgeBaseM{ss} = phiM;  edgeBaseP{ss} = phiP;
end