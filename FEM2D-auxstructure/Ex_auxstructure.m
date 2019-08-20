clc;clear;close all;

msh = 'triangulation';
if strcmp(msh,'triangulation') % triangulation
    [node,elem] = squaremesh([0 1 0 1],0.5);
    showmesh(node,elem);
    findnode(node); findelem(node,elem); findedge(node,elem);
else  % polygonal meshes
    load('meshex1.mat');
    showmesh(node,elem);
    findnode(node); findelem(node,elem); findedge(node,elem);
end

aux = auxstructure(node,elem)

