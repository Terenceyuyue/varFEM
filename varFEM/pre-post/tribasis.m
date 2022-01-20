function Dtri =  tribasis(i,j,k,Dlambda,lambda)
%tribasis computes second derivatives of lambda_i*lambda_j*lambda_k
% note: Dtri(:,:,1) = D11, Dtri(:,:,2) = D22, Dtri(:,:,3) = D12
%
% Copyright (C) Terence Yu.
    
    NT = size(Dlambda,1); nQuad = size(lambda,1);
    Dtri = zeros(NT,nQuad,3); % [11,22,12]
    ss  = [1 2 1]; tt = [1 2 2]; 
    for m = 1:3  % loop for [11,22,12]
        s = ss(m); t = tt(m);
        b1 = (Dlambda(:,s,i).*Dlambda(:,t,j)+Dlambda(:,t,i).*Dlambda(:,s,j))*lambda(:,k)';
        b2 = (Dlambda(:,s,j).*Dlambda(:,t,k)+Dlambda(:,t,j).*Dlambda(:,s,k))*lambda(:,i)';
        b3 = (Dlambda(:,s,i).*Dlambda(:,t,k)+Dlambda(:,t,i).*Dlambda(:,s,k))*lambda(:,j)';
        Dtri(1:NT,:,m) = b1 + b2 + b3; 
    end
end