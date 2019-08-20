clc;clear;close all;

load('meshex1.mat');

% showmesh
figure,
showmesh(node,elem);
findnode(node);
findelem(node,elem);

% showsolution for u = sin(x)*cos(y)
figure,
x = node(:,1); y = node(:,2); u = sin(x).*cos(y);
showsolution(node,elem,u);

% show the solution by using showmesh
figure,
data = [node,u];
showmesh(data,elem);
