function I = integral3d(Th,fh,Vh,quadOrder)
%  int_\Omega f(x,y,z) dxdy, where f(x,y,z) is apparoximated by the FEM
%  function:
%    f(x,y,z) \approx fh(x,y,z) = f1*Phi1(x,y,z) + f2*Phi2(x,y,z) + ....
%  Then,
%    I = f1 * \int_\Omega(Phi1) + f2 * \int_\Omega(Phi2) + ....
%    
%  fh: - a function handle
%      - an FEM function 
%      - Coef matrix
%

I = assem3d(Th,fh,[],[],Vh,quadOrder);