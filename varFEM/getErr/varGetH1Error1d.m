function ErrH1 = varGetH1Error1d(Th,Due,uh,Vh,quadOrder)

%Due = pde.Du;  % Du = u_x

% ErrH1
uxe = @(p) Due(p);  
uxec = interp1dMat(uxe,'.val',Th,Vh,quadOrder);
uhxec = interp1dMat(uh,'.dx',Th,Vh,quadOrder);

cc = (uxec-uhxec).^2 ;
ErrH1 = integral1d(Th,cc,Vh,quadOrder);
ErrH1 = sqrt(ErrH1);