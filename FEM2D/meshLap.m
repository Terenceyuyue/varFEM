function [node,elem,eD,elemN] = meshLap()

a1 = 0; b1 = 1; a2 = 0; b2 = 1; Nx = 5; Ny = 5;
square = [a1 b1 a2 b2]; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh(square,h1,h2);
showmesh(node,elem); findnode(node); findelem(node,elem);

eps = 1e-4; 
% --------- Set Dirichlet boundaries ---------------
idL = find(abs(node(:,1)-a1)<=eps); % left
idD = find(abs(node(:,2)-a2)<=eps); % down
idU = find(abs(node(:,2)-b2)<=eps); % upper
id = [idL;idD;idU]; 
eD = unique(id); % node id of Dirichlet boundaries

% --------- Set Neumann boundaries ---------------
idL = find(abs(node(:,1)-b1)<=eps); % right
yn = node(idL,2); 
[~,m] = sort(yn); % sort in ascending order
id = idL(m);
% Neumann elements
nel = length(id)-1; % number of Nuemann boundaries
elemN = zeros(nel,2); 
% id of starting and ending points
elemN(:,1) = id(1:end-1); elemN(:,2) = id(2:end); 
