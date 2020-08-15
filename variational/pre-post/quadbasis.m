function Dquad = quadbasis(i,j,Dphi,lambda)
%quadbasis computes second derivatives of lambda_i*lambda_j
% note: Dquad(:,:,1) = D11, Dquad(:,:,2) = D22, Dquad(:,:,3) = D12
%
% Copyright (C) Terence Yu.

    NT = size(Dphi,1); nQuad = size(lambda,1);
    Dquad = zeros(NT,nQuad,3); % [11,22,12]
    Dquad(1:NT,:,1) = repmat(2*Dphi(:,1,i).*Dphi(:,1,j),1,nQuad);
    Dquad(1:NT,:,2) = repmat(2*Dphi(:,2,i).*Dphi(:,2,j),1,nQuad);
    Dquad(1:NT,:,3) = repmat((Dphi(:,1,i).*Dphi(:,2,j) + Dphi(:,2,i).*Dphi(:,1,j)), 1,nQuad);
end