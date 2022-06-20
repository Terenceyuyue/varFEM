clc;clear;close all

BdBox = [0 1 0 1]; h0 = 0.1;
[x,y] = meshgrid(BdBox(1):h0:BdBox(2),BdBox(3):h0:BdBox(4));  
p = [x(:),y(:)];  

t = my_delaunay(p);
showmesh(p,t)
