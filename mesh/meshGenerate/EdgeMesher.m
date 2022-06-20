function [node,elem] = EdgeMesher(P1,P2,options)
%EdgeMesher Generate triangular meshes according to the line segments for
% the boundary.
% The line segments must be in the conventional directions.
%
%    [node,elem] = EdgeMesher(P1,P2)
%    P1 = [x,y] is the starting points of edges
%    P2 = [x,y] is the ending points of edges
%    
%    options: plot is true or false
%             refineTimes is the number of refinements
%


if nargin==2
    options.plot = false; 
end

%% Decomposed geometry matrix
N = size(P1,1);
g = zeros(7,N);
g(1,:) = 2;  % 2 for line 
g(2,:) = P1(:,1); % x: starting points
g(3,:) = P2(:,1); % x: ending points
g(4,:) = P1(:,2); % y: starting points
g(5,:) = P2(:,2); % y: ending points
g(6,:) = 1; % label of subdomain on the left
g(7,:) = 0; % label of subdomain on the right

%% Triangulation
if ~isfield(options,'hmax')
    options.hmax = 1;
end
[p,e,t] = initmesh(g,'hmax',options.hmax); % initial mesh
if ~isfield(options,'refineTimes')
    options.refineTimes = 1; 
end
for i = 1:options.refineTimes
    [p,e,t] = refinemesh(g,p,e,t);
end
node = p';
elem = t(1:3,:)';


%% Plot mesh
if ~isfield(options,'plot')
    options.plot = false;
end
if options.plot
    figure,
    showmesh(node,elem); 
end