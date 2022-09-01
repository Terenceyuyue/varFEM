function ErrL2 = varGetL2Error1d(Th,ue,uh,Vh,quadOrder)

%ue = pde.uexact

% ErrL2
uec = interp1dMat(ue,'.val',Th,Vh,quadOrder);
uhc = interp1dMat(uh,'.val',Th,Vh,quadOrder);
cc = (uec-uhc).^2;
ErrL2 = integral1d(Th,cc,Vh,quadOrder);
ErrL2 = sqrt(ErrL2);