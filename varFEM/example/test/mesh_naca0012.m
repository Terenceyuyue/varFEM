function [node,elem] = mesh_naca0012()

%% Circle
N = 50; r = 5;
theta = linspace(0,2*pi,N+1)';
x = r*cos(theta);  
y = r*sin(theta);  
PC1 = [x(1:end-1), y(1:end-1)];
PC2 = [x(2:end), y(2:end)];
PC2(end,:) = [x(1), y(1)];

%% naca0012
naca12 = @(x) 0.17735*sqrt(x) - 0.075597*x - 0.212836*(x.^2) ...
    + 0.17363*(x.^3) - 0.06254*(x.^4);
N = 90;
x = linspace(0,1,N+1)';
% lower
x1 = x(end:-1:1);
y1 = -naca12(x1);
% upper
x2 = x;
y2 = naca12(x2);
% combine
x = [x1; x2(2:end-1)];
y = [y1; y2(2:end-1)];
Pnaca1 = [x(1:end-1), y(1:end-1)];
Pnaca2 = [x(2:end), y(2:end)];
Pnaca2(end,:) = [x(1), y(1)];

%% air region
P1 = [PC1; Pnaca1];
P2 = [PC2; Pnaca2];
options.refineTimes = 0;
[node,elem] = EdgeMesher(P1,P2,options);

% figure, showmesh(node,elem);
% x = linspace(0,1,1000)';
% hold on
% plot(x,naca12(x),'r-','linewidth',1);
% plot(x,-naca12(x),'r-','linewidth',1);
% hold off
% axis([-0.5 1.5 -0.5 0.5]);