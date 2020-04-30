function bdStruct= setboundary1D(node,elem1D,varargin)
% varargin: string for Neumann boundary

% edge number
bdEdge = find(accumarray(elem1D(:),1) == 1);
nE = size(bdEdge,1); 

% initial as Dirichlet (true for Dirichlet, false for Neumann)
bdFlag = true(nE,1);  nodeEdge = node(bdEdge,:);
if size(node,2)==1, x = nodeEdge; end
if size(node,2)==2, x = nodeEdge(:,1); y = nodeEdge(:,2); end
nvar = length(varargin); % 1 * size(varargin,2)
% note that length(varargin) = 1 for bdNeumann = [] or ''
if (nargin==2) || (~isempty(varargin{1})) 
    for i = 1:nvar
        bdNeumann = varargin{1};
        id = eval(bdNeumann);
        bdFlag(id) = false;
    end
end
bdStruct.Dirichlet = bdEdge(bdFlag); % Dirichlet boundary nodes
bdStruct.Neumann = bdEdge(~bdFlag); % Dirichlet boundary nodes