function bdStruct= setboundary(node,elem,varargin)
% varargin: string for Neumann boundary

% ------- totalEdge ---------
if size(elem,2) == 4
    elem = mat2cell(elem,ones(size(elem,1),1), 4); % rectangle element
end
if iscell(elem)
    shiftfun = @(verts) [verts(2:end),verts(1)];  % or shiftfun = @(verts) circshift(verts,-1);
    T1 = cellfun(shiftfun, elem, 'UniformOutput', false);
    v0 = horzcat(elem{:})'; % the starting points of edges
    v1 = horzcat(T1{:})'; % the ending points of edges
    allEdge = [v0,v1];
else % Triangulation
    allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
end
totalEdge = sort(allEdge,2);

% --------  counterclockwise bdEdge --------
[~,~,s] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
[~, i1, ~] = unique(totalEdge,'rows');
bdEdge = allEdge(i1(s==1),:);

% --------- set boundary --------
nE = size(bdEdge,1);
% initial as Dirichlet (true for Dirichlet, false for Neumann)
bdFlag = true(nE,1);
nodebdEdge = (node(bdEdge(:,1),:) + node(bdEdge(:,2),:))/2;
x = nodebdEdge(:,1); y = nodebdEdge(:,2);
nvar = length(varargin); % 1 * size(varargin,2)
% note that length(varargin) = 1 for bdNeumann = [] or ''
if (nargin==2) || (~isempty(varargin{1})) 
    for i = 1:nvar 
        bdNeumann = varargin{i};
        id = eval(bdNeumann);
        bdFlag(id) = false;
    end
end

bdStruct.elemD = bdEdge(bdFlag,:); % Dirichlet boundary edges
bdStruct.elemN = bdEdge(~bdFlag,:); % Neumann boundary edges
bdStruct.eD = unique(bdEdge(bdFlag,:)); % Dirichlet boundary nodes
bdIndex = find(s==1);      % indices of all boundary edges
bdStruct.bdIndex = bdIndex; 
bdStruct.bdIndexD = bdIndex(bdFlag); % indices of Dirichelt boundary edges
bdStruct.bdIndexN = bdIndex(~bdFlag); % indices of Neumann boundary edges