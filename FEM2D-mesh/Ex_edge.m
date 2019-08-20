clc;clear;close all;

% % -------------- edge (triangulation) ----------
% [node,elem] = squaremesh([0 1 0 1],0.5);
% figure, % showmesh
% showmesh(node,elem); 
% findelem(node,elem); findnode(node);
% 
% figure, % boundary edges
% showmesh(node,elem); 
% bdInd = 1;
% findedge(node,elem,bdInd);
% figure, % all edges
% showmesh(node,elem); 
% findedge(node,elem);

% -------------- edge (polygonal meshes) ----------
load('meshex1.mat');
figure, % showmesh
showmesh(node,elem);
findelem(node,elem); findnode(node);

figure, % boundary edges
showmesh(node,elem);
bdInd = 1;
findedge(node,elem,bdInd);
figure, % all edges
showmesh(node,elem);
findedge(node,elem);
