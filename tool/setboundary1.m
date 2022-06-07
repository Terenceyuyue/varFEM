function bdStruct= setboundary1(node,elem1,varargin)
% varargin: string for Neumann boundary
% elem1 may correspond to 2-D line segments

% edge number
bdEdge = find(accumarray(elem1(:),1) == 1);
nE = size(bdEdge,1); 

% initialize as Dirichlet (true for Dirichlet, false for Neumann)
bdNodeIdxD = true(nE,1);  nodeEdge = node(bdEdge,:);
if size(node,2)==1, x = nodeEdge; end %#ok<NASGU> 
if size(node,2)==2, x = nodeEdge(:,1); y = nodeEdge(:,2); end %#ok<NASGU> 
nvar = length(varargin); % 1 * size(varargin,2)
% note that length(varargin) = 1 for bdNeumann = [] or ''
if (nargin==2) || (~isempty(varargin{1})) 
    for i = 1:nvar
        bdNeumann = varargin{i};
        if mycontains(bdNeumann,'==')    % transform 'x==1' to 'abs(x-1)<1e-4'
            bd = bdNeumann;
            bd = strrep(bd,'==','-');
            if mycontains(bd,'|')
                bd = strrep(bd, '|', ')<1e-4 | abs(');
            end
            bdNeumann = ['abs(',   bd,   ')<1e-4'];
        end
        id = eval(bdNeumann);
        bdNodeIdxD(id) = false;
    end
end
bdStruct.Dirichlet = bdEdge(bdNodeIdxD); % Dirichlet boundary nodes
bdStruct.Neumann = bdEdge(~bdNodeIdxD); % Neumann boundary nodes