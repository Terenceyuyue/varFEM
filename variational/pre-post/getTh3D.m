function Th = getTh3D(node,elem,bdNeumann)
%% GETTH3D gets 3D mesh information

if nargin==2, bdNeumann = []; end

% node, elem
elem = fixorder3(node,elem);  
% % WARNING: the output of fixorder3.m can not be used in uniformrefine3.m
Th.node = node; Th.elem = elem;

% bdStruct
Th.bdStruct = setboundary3(node,elem,bdNeumann);

% edge, elem2edge
allEdge = [elem(:,[1,2]); elem(:,[1,3]); elem(:,[1,4]);
           elem(:,[2,3]); elem(:,[2,4]); elem(:,[3,4])];
totalEdge = sort(allEdge,2);
[edge, ~, totalJ] = unique(totalEdge,'rows');
Th.edge = edge;
NT = size(elem,1);
Th.elem2edge = reshape(totalJ, NT, 6);

% numbers
Th.N = size(node,1); Th.NT = NT; 
Th.NE = size(Th.edge,1);
