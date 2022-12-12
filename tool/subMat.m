function Aij = subMat(A,I,J)
%Get (I,J) entries of A, where I and I are row and column subscripts with the
%same size.

Aij = A(sub2ind(size(A),I,J));
