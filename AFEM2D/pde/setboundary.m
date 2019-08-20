function bdFlag = setboundary(node,a1,b1,a2,b2)

eps = 1e-4; 
% --------- Set Dirichlet boundaries ---------------
idL = find(abs(node(:,1)-a1)<=eps); % left
idR = find(abs(node(:,1)-b1)<=eps); % right
idD = find(abs(node(:,2)-a2)<=eps); % down
idU = find(abs(node(:,2)-b2)<=eps); % upper
id = [idL;idR;idD;idU]; 
eD = unique(id); % node id of Dirichlet boundaries

% --------- Set Neumann boundaries ---------------
elemN = [];

bdFlag.eD = eD; bdFlag.elemN = elemN;