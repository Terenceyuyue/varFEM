clc;clear;close all

[node,elem] = squaremesh([0,1,0,1],0.5);
showmesh(node,elem); 
findnode(node); findelem(node,elem); findedge(node,elem);

bdNeumann = 'abs(x-0)<1e-4';
bdStruct = setboundary(node,elem,bdNeumann);