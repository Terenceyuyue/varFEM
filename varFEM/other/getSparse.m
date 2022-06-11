function [ii,jj] = getSparse(elem2dofv, elem2dofu)
%%getSparse returns the sparse indices of Aij in block matrix
%

NT = size(elem2dofv,1);
Ndofv = size(elem2dofv,2); Ndofu = size(elem2dofu,2);

nnz = NT*Ndofv*Ndofu;
ii = zeros(nnz,1); jj = zeros(nnz,1); id = 0;
for i = 1:Ndofv
    for j = 1:Ndofu
        ii(id+1:id+NT) = elem2dofv(:,i);   
        jj(id+1:id+NT) = elem2dofu(:,j);   
        id = id + NT;
    end
end