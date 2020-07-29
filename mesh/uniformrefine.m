function [node,elem] = uniformrefine(node,elem)
%uniformrefine uniformly refines a 2-D triangulation

% Copyright (C) Long Chen, modified by Terence Yu. 

%% auxiliary mesh data
auxT = auxstructure(node,elem);
edge = auxT.edge; elem2edge = auxT.elem2edge;
N = size(node,1); NT = size(elem,1); NE = size(edge,1); 

%% Add new nodes: middle points of all edges
node(N+1:N+NE,:) = (node(edge(:,1),:)+node(edge(:,2),:))/2;

%% Refine each triangle into four triangles as follows
% 3
% | \
% 5- 4
% |\ |\
% 1- 6- 2
t = 1:NT; p = zeros(NT,6);
p(:,1:3) = elem;
p(:,4:6) = elem2edge + N;
elem(t,:) = [p(t,1), p(t,6), p(t,5)];
elem(NT+1:2*NT,:) = [p(t,6), p(t,2), p(t,4)];
elem(2*NT+1:3*NT,:) = [p(t,5), p(t,4), p(t,3)];
elem(3*NT+1:4*NT,:) = [p(t,4), p(t,5), p(t,6)];
