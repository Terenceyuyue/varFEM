function Th = getTh(node,elem,bdNeumann)

if nargin==2, bdNeumann = []; end
% node, elem
Th.node = node; Th.elem = elem; 
% bdStruct
bdStruct = setboundary(node,elem, bdNeumann);
Th.bdStruct = bdStruct;

% aux,auxT
aux = auxgeometry(node,elem);  auxT = auxstructure(node,elem);
Th.aux = aux; Th.auxT = auxT;
edge = auxT.edge;

% numbers
Th.N = size(node,1); Th.NT = size(elem,1); Th.NE = size(edge,1);