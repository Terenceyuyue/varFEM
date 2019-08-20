%FEM1D solves the 1-D partial diffential equation
%  
%  au''+bu+cu = f,  x \in (x0,xL);
%     g_N = du(x0), g_D = u(xL);
%  or g_D = u(x0),  g_N = du(xL);
%
% where, g_N for Neumann boundary conditions and 
%        g_D for Dirichlet boundary conditions.
%
% Copyright (C) Terence YUE Yu. 

clc;clear; close all;
tic;
% ------ Mesh -----------
a = 0; b = 1;
nel = 10;  N = nel+1; % numbers of elements and nodes
node = linspace(a,b,nel+1)';
elem = zeros(nel,2); elem(:,1) = 1:N-1; elem(:,2) = 2:N;

% ------ PDE ---------
acoef = 1;  bcoef = 1;  ccoef = 1;
syms x;
c1 = 0.5/exp(1); c2 = -0.5*(1+1/exp(1));
u = c1*exp(2*x)+c2*exp(x)+1/2;
% exact solution
uexact = eval(['@(x)',vectorize(u)]);  % transform to anonymous function

du = diff(u);
du = eval(['@(x)',vectorize(du)]);

% rhs
f = acoef*diff(u,2)+bcoef*diff(u,1)+ccoef*u;
f = eval(['@(x)',vectorize(f)]); 

% boundary conditions
Neumann = 1; Dirichlet = N;
g_N = du(node(Neumann)); 
g_D = uexact(node(Dirichlet));

% --------- Assembly ----------
kk = zeros(N,N);  ff = zeros(N,1);
for iel = 1:nel
    % local --> global
    index = elem(iel,:);
    nl = index(1); nr = index(2);
    xl = node(nl); xr = node(nr); he = xr-xl;
    
    % element matrix and vector
    a1 = -(acoef/he);  a2 = bcoef/2;  a3 = ccoef*he/6;
    ke = a1*[1 -1;-1 1]+a2*[-1 1;-1 1]+a3*[2 1;1 2];
    fe = f((xl+xr)/2)*[he/2;he/2];
    
    % assemble
    kk(index,index) = kk(index,index)+ke;
    ff(index) = ff(index)+fe;
end
kk = sparse(kk);

% ------- Neumann boundary conditions -----------
if ~isempty(Neumann)
    bn = [-1;zeros(N-2,1);1]; % -1: left norm vector
    bn(Neumann) = bn(Neumann)*g_N;
    ff = ff + (-acoef)*bn;
end

% ------- Dirichlet boundary conditions ----------
if ~isempty(Dirichlet)
    isBdNode = false(N,1); isBdNode(Dirichlet) = true;
    bdNode = find(isBdNode); freeNode = find(~isBdNode);
    u = zeros(N,1); u(bdNode) = g_D;
    ff = ff - kk*u;
end

% -------- error analysis -----------
uexact = uexact(node);
u(freeNode) = kk(freeNode,freeNode)\ff(freeNode);
figure,plot(node,u,'r-',node,uexact,'k--','linewidth',1)
xlabel('x');ylabel('u')
legend('Numerical solution','Exact solution')
Err = u-uexact;
figure, plot(node,Err,'linewidth',1)
legend('Absolute error')
toc