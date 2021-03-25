function Dquad = quadbasis(i,j,Dlambda,lambda)
%quadbasis computes second derivatives of lambda_i*lambda_j
% note: Dquad(:,:,1) = D11, Dquad(:,:,2) = D22, Dquad(:,:,3) = D12
%
% Copyright (C) Terence Yu.

    NT = size(Dlambda,1); nQuad = size(lambda,1);
    Dquad = zeros(NT,nQuad,3); % [11,22,12]
    Dquad(1:NT,:,1) = repmat(2*Dlambda(:,1,i).*Dlambda(:,1,j),1,nQuad);
    Dquad(1:NT,:,2) = repmat(2*Dlambda(:,2,i).*Dlambda(:,2,j),1,nQuad);
    Dquad(1:NT,:,3) = repmat((Dlambda(:,1,i).*Dlambda(:,2,j) + Dlambda(:,2,i).*Dlambda(:,1,j)), 1,nQuad);
end