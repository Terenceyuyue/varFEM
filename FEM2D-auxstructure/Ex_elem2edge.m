clc;clear;close all;

% % -------- elem2edge (triangulation) --------
% [node,elem] = squaremesh([0 1 0 1],0.5);
% showmesh(node,elem); findelem(node,elem); findnode(node);
% findedge(node,elem);

% -------- elem2edge (polygonal meshes) --------
load('meshex1.mat');
showmesh(node,elem);
findnode(node); findelem(node,elem); findedge(node,elem);

NT = size(elem,1);

if iscell(elem)
    % totalEdge
    shiftfun = @(verts) [verts(2:end),verts(1)];  % or shiftfun = @(verts) circshift(verts,-1);
    T1 = cellfun(shiftfun, elem, 'UniformOutput', false);
    v0 = horzcat(elem{:})'; % the starting points of edges
    v1 = horzcat(T1{:})'; % the ending points of edges
    totalEdge = sort([v0,v1],2);
    
    % elem2edge
    [~, ~, totalJ] = unique(totalEdge,'rows');
    cellLen = cellfun('length',elem); % length of each cell
    elem2edge = mat2cell(totalJ,cellLen,1);
    elem2edge = cellfun(@transpose, elem2edge, 'UniformOutput', false);
    
else % Triangulation
    totalEdge = sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2);
    [~, ~, totalJ] = unique(totalEdge,'rows');    
    elem2edge = reshape(totalJ,NT,3);
end

% -------- edge, bdEdge --------
[i,j,s] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
edge = [j,i];
bdEdge = edge(s==1,:);


