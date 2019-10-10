function D = Dquad(i,j,iel,Dphi,lambda)

Deriv = @(s,t) Dphi(iel,s,i).*Dphi(iel,t,j) ...
                + Dphi(iel,t,i).*Dphi(iel,s,j);

D = [Deriv(1,1),Deriv(2,2),Deriv(1,2)];
nI = size(lambda,1);
D = repmat(D,nI,1);
