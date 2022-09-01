function ErrL2 = varGetL2Error(Th,ue,uh,Vh,quadOrder)

%ue = pde.uexact

% ErrL2
uec = interp2dMat(ue,'.val',Th,Vh,quadOrder);
uhc = interp2dMat(uh,'.val',Th,Vh,quadOrder);
cc = (uec-uhc).^2;
ErrL2 = integral2d(Th,cc,Vh,quadOrder);
ErrL2 = sqrt(ErrL2);