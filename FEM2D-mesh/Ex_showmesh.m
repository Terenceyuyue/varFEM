clc;clear;close all;

% -------- showmesh - meshex1.mat --------

load('meshex1.mat');

% max_n_vertices = max(cellfun(@length, elem));
% padding_func = @(vertex_ind) [vertex_ind,...
%       NaN(1,max_n_vertices-length(vertex_ind))];  % function to pad the vacancies (∫·œÚ∆¥Ω”)
% tpad = cellfun(padding_func, elem, 'UniformOutput', false);
% tpad = vertcat(tpad{:});
% h = patch('Faces', tpad,'Vertices', node);
% set(h,'facecolor',[0.5 0.9 0.45],'edgecolor','k');
% axis equal; axis tight;

showmesh(node,elem);
findnode(node);
findelem(node,elem);