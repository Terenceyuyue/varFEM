clc;clear;close all

BdBox = [0 1 0 1]; h0 = 0.25;
[x,y] = meshgrid(BdBox(1):h0:BdBox(2),BdBox(3):h0:BdBox(4));  
node = [x(:),y(:)];  
id = randperm(size(node,1));
node = node(id,:);
elem = my_delaunay(node);

figure,
subplot(1,2,1), showmesh(node,elem);
subplot(1,2,2), 
signSwap = -1;
[node,elem] = edgeswap(node,elem,signSwap);
showmesh(node,elem)
findnode(node)
% findelem(node,elem)
% findedge(node,elem)
