function [kk,ff] = Get_kk_ff(node,elem,pde)

N = size(node,1);  NT = size(elem,1); Ndof = 3; 

% area of all triangules
x1 = node(elem(:,1),1); y1 = node(elem(:,1),2);
x2 = node(elem(:,2),1); y2 = node(elem(:,2),2);
x3 = node(elem(:,3),1); y3 = node(elem(:,3),2);
A = 0.5*((x2.*y3-x3.*y2)-(x1.*y3-x3.*y1)+(x1.*y2-x2.*y1)); 

% All element matrices
k11 = (x3-x2).^2+(y2-y3).^2; 
k12 = (x3-x2).*(x1-x3)+(y2-y3).*(y3-y1);
k13 = (x3-x2).*(x2-x1)+(y2-y3).*(y1-y2);
k21 = k12; 
k22 = (x1-x3).^2+(y3-y1).^2;
k23 = (x1-x3).*(x2-x1)+(y3-y1).*(y1-y2);
k31 = k13;
k32 = k23;
k33 = (x2-x1).^2+(y1-y2).^2;
K = [k11,k12,k13,k21,k22,k23,k31,k32,k33];
K = K./(4*repmat(A,1,Ndof^2));

% All vectors
f = pde.f;
xc = 1/3*(x1+x2+x3); yc = 1/3*(y1+y2+y3); pc = [xc,yc];
F1 = f(pc).*A./3; F2 = F1; F3 = F1;
F = [F1,F2,F3];

% local --> global
nnz = NT*Ndof^2;
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = K(:);
id = 0; 
for i = 1:Ndof
    for j = 1:Ndof
        ii(id+1:id+NT) = elem(:,i);   % zi
        jj(id+1:id+NT) = elem(:,j);   % zj
        id = id + NT; 
    end
end

% stiff matrix and load vector
kk = sparse(ii,jj,ss,N,N);
ff = accumarray(elem(:), F(:),[N 1]);
end

