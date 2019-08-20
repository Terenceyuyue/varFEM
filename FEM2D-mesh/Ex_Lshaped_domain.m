clc;clear;close all;
% -------- L-shaped domain --------
node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1; 0,-1; 0,0];  % coordinates
elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];     % connectivity
showmesh(node,elem);
grid off; box on;
findelem(node,elem);  % plot indices of all triangles
findnode(node);       % plot indices of all vertices