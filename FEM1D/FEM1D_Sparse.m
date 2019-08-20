%FEM1D_Sparse solves the 1-D partial diffential equation 
%  
%  au''+bu+cu = f,  x \in (x0,xL);
%     g_N = du(x0), g_D = u(xL);
%  or g_D = u(x0),  g_N = du(xL);
%
% where, g_N for Neumann boundary conditions and 
%        g_D for Dirichlet boundary conditions.
% 
% The sparse assembly is frequently used in the FEM programming.
%
% Copyright (C) Terence YUE Yu. 

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

% --------------------- Assemble the matrix ----------------
% All element matrices 
h = diff(node);
k11 = -acoef./h+bcoef/2*(-1)+ccoef*h./6*2;      
k12 = -acoef./h*(-1)+bcoef/2+ccoef*h./6;
k21 = -acoef./h*(-1)+bcoef/2*(-1)+ccoef*h./6;   
k22 = -acoef./h+bcoef/2+ccoef*h./6*2;
K = [k11,k12,k21,k22]; % stored in rows

% local --> global
Ndof = 2; nnz = nel*Ndof^2;
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = K(:);
id = 0; 
for i = 1:Ndof
    for j = 1:Ndof
        ii(id+1:id+nel) = elem(:,i);   % zi
        jj(id+1:id+nel) = elem(:,j);   % zj
        id = id + nel; 
    end
end

% stiff matrix
kk = sparse(ii,jj,ss,N,N);

% --------- Assemble the vector ------------
x1 = node(1:N-1); x2 = node(2:N);
xc = (x1+x2)./2; 
F1 = f(xc).*h./2; F2 = F1; F = [F1,F2];
ff = accumarray(elem(:), F(:),[N 1]);

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