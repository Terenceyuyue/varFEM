function Dtri = Dtriple(i,j,k,iel,Dphi,lambda)

Deriv = @(s,t) (Dphi(iel,s,i).*Dphi(iel,t,j)+Dphi(iel,t,i).*Dphi(iel,s,j))*lambda(:,k) ...
            + (Dphi(iel,s,j).*Dphi(iel,t,k)+Dphi(iel,t,j).*Dphi(iel,s,k))*lambda(:,i) ...
            + (Dphi(iel,s,i).*Dphi(iel,t,k)+Dphi(iel,t,i).*Dphi(iel,s,k))*lambda(:,j);

% nI = size(lambda,1); % number of integration points        
% Dtri = zeros(nI,3); % (11,22,12)
Dtri = [Deriv(1,1), Deriv(2,2), Deriv(1,2)];