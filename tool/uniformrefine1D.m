function [node,elem1D] = uniformrefine1D(node,elem1D)

N = size(node,1); nel = size(elem1D,1); ndof = 2;
% add new nodes
node(N+1:N+nel) = (node(elem1D(:,1))+node(elem1D(:,2)))/2;

% add new elements
% 1 --- 3 --- 2
t = 1:nel; p = zeros(nel,2*ndof);
p(:,1:2) = elem1D;  p(:,3) = (1:nel)' + N;
elem1D(t,:) = [p(t,1), p(t,3)];
elem1D(nel+1:2*nel,:) = [p(t,3), p(t,2)];
