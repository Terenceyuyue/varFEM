function ErrL2 = varGetL2Error3d(Th,ue,uh,Vh,quadOrder)

%ue = pde.uexact

% ErrL2
uec = interp3dMat(ue,'.val',Th,Vh,quadOrder);
uhc = interp3dMat(uh,'.val',Th,Vh,quadOrder);
cc = (uec-uhc).^2;
ErrL2 = integral3d(Th,cc,Vh,quadOrder);
ErrL2 = sqrt(ErrL2);