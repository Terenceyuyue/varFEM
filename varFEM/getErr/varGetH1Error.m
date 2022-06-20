function ErrH1 = varGetH1Error(Th,Due,uh,Vh,quadOrder)

%Due = pde.Du;

% ErrH1
uxe = @(p) Due(p)*[1;0];  uye = @(p) Due(p)*[0;1];
uxec = interp2dMat(uxe,'.val',Th,Vh,quadOrder);
uhxec = interp2dMat(uh,'.dx',Th,Vh,quadOrder);
uyec = interp2dMat(uye,'.val',Th,Vh,quadOrder);
uhyec = interp2dMat(uh,'.dy',Th,Vh,quadOrder);
cc = (uxec-uhxec).^2 + (uyec-uhyec).^2;
ErrH1 = integral2d(Th,cc,Vh,quadOrder);
ErrH1 = sqrt(ErrH1);