function Th = getTh3D(node,elem,bdNeumann)
%getTh3D gets 3D mesh information
%
% Copyright (C) Terence Yu.

if nargin==2, bdNeumann = []; end

% node, elem
elem = fixorder3(node,elem);  % with reverse normal vector 
elem(:,[2,3]) = elem(:,[3,2]);

% % WARNING: the output of fixorder3.m can not be used in uniformrefine3.m
Th.node = node; Th.elem = elem;

% bdStruct
bdStruct = setboundary3(node,elem,bdNeumann);

% edge, elem2edge, face2edge
Th.edge = bdStruct.edge;
Th.elem2edge = bdStruct.elem2edge;
Th.face2edge = bdStruct.face2edge;

bdStruct = rmfield(bdStruct,{'edge','elem2edge','face2edge'});
Th.bdStruct = bdStruct;

% numbers
Th.N = size(node,1); Th.NT = size(elem,1); 
Th.NE = size(Th.edge,1);
