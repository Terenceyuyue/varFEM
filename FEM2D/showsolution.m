function showsolution(node,elem,u)
%Showsolution displays the solution corresponding to a mesh given by [node,elem] in 2-D.
%
% Copyright (C) Terence Yue Yu.

data = [node,u];
patch('Faces', elem,...
    'Vertices', data,...
    'FaceColor', 'interp',...
    'CData', u / max(abs(u)) );
axis('square');
sh = 0.05;
xlim([min(node(:,1)) - sh, max(node(:,1)) + sh])
ylim([min(node(:,2)) - sh, max(node(:,2)) + sh])
zlim([min(u) - sh, max(u) + sh])
xlabel('x'); ylabel('y'); % zlabel('u');

view(3); grid on; % view(150,30);  
