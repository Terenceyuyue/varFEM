clc;clear; close all;

tic;
% -------------- Mesh -----------
a = 0; b = 1;
nel = 10;  N = nel+1; % numbers of elements and nodes
node = linspace(a,b,nel+1)';
elem = zeros(nel,2); elem(:,1) = 1:N-1; elem(:,2) = 2:N;

% --------------- PDE ------------
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

% ----- All element matrices --------
h = diff(node);
k11 = -acoef./h+bcoef/2*(-1)+ccoef*h./6*2;      k12 = -acoef./h*(-1)+bcoef/2+ccoef*h./6;
k21 = -acoef./h*(-1)+bcoef/2*(-1)+ccoef*h./6;   k22 = -acoef./h+bcoef/2+ccoef*h./6*2;

% -------- local --> global ---------
z1 = elem(:,1); z2 = elem(:,2);

% ----------- Assemble the matrix ----------
% upper triangular
kk = sparse(z1,z2,k12,N,N); 
% lower triangular
kk = kk+sparse(z2,z1,k21,N,N); 
% diagonal
kk = kk+sparse(z1,z1,k11,N,N);
kk = kk+sparse(z2,z2,k22,N,N);

% --------- Assemble the vector ------------
x1 = node(1:N-1); x2 = node(2:N); xc = (x1+x2)./2;
fi = f(xc).*h./2;
F1 = fi; F2 = fi;
ff = sparse(z1,1,F1,N,1);
ff = ff+sparse(z2,1,F2,N,1);

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

