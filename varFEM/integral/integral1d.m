function I = integral1d(Th,fh,Vh,quadOrder)
%  int_\Omega f(x) dx, where f(x) is apparoximated by the FEM
%  function:
%    f(x) \approx fh(x) = f1*Phi1(x) + f2*Phi2(x) + ....
%  Then,
%    I = f1 * \int_\Omega(Phi1) + f2 * \int_\Omega(Phi2) + ....
%    
%  fh: - a function handle
%      - an FEM function 
%      - Coef matrix
%

I = assem1d(Th,fh,[],[],Vh,quadOrder);