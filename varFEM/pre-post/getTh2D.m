function Th = getTh2D(node,elem,bdNeumann)
%% GETTH2D gets 2D mesh information
%
% Copyright (C) Terence Yu.

if nargin==2, bdNeumann = []; end

% node, elemD
Th.node = node; Th.elem = elem;

% bdStruct
Th.bdStruct = setboundary(node,elem,bdNeumann);

% edge, elem2edge
allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
totalEdge = sort(allEdge,2);
[edge, ~, totalJ] = unique(totalEdge,'rows');
NT = size(elem,1);
elem2edge = reshape(totalJ,NT,3);
Th.edge = edge; 
Th.elem2edge = elem2edge;

% numbers
Th.N = size(node,1); Th.NT = NT; Th.NE = size(edge,1);

end
