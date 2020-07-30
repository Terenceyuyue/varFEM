function Th = getTh1D(node,elem,bdNeumann)
%% GETTH1D gets 1D mesh information


% node, elem1D
Th.node = node; Th.elem = elem; 

% bdStruct
Th.bdStruct = setboundary1(node,elem,bdNeumann);

% numbers
Th.N = size(node,1); Th.nel = size(elem,1);

end