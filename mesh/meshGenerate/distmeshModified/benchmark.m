clc;clear;close all

num = 3;
if num==1
    disp('Example 1: Unit Circle');
    fd = @(p) sqrt(sum(p.^2,2))-1;
    h0 = 0.2;
    [node,elem]=distmesh(fd,h0,[-1 1 -1 1]);
    figure, showmesh(node,elem)
end

if num==2
    disp('Example 2: Square with hole');
    fd = @(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5));
    h0 = 0.2;
    BdBox = [-1 1 -1 1];
    pfix = [-1,-1;-1,1;1,-1;1,1];
    [node,elem]=distmesh(fd,h0,BdBox,pfix);
    figure, showmesh(node,elem)
end

if num==3
    disp('Example 4: Ellipse');
    fd = @(p) p(:,1).^2/2^2+p(:,2).^2/1^2-1;
    BdBox = [-2 2 -1 1];
    h0 = 0.1;    
    [node,elem] = distmesh(fd,h0,BdBox);
    figure, showmesh(node,elem)
end

