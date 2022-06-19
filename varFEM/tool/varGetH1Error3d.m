function ErrH1 = varGetH1Error3d(Th,Due,uh,Vh,quadOrder)

%Due = pde.Du;

% ErrH1
uxe = @(p) Due(p)*[1;0;0];  
uye = @(p) Due(p)*[0;1;0]; 
uze = @(p) Due(p)*[0;0;1];

uxec = interp3dMat(uxe,'.val',Th,Vh,quadOrder);
uyec = interp3dMat(uye,'.val',Th,Vh,quadOrder);
uzec = interp3dMat(uze,'.val',Th,Vh,quadOrder);

uhxec = interp3dMat(uh,'.dx',Th,Vh,quadOrder);
uhyec = interp3dMat(uh,'.dy',Th,Vh,quadOrder);
uhzec = interp3dMat(uh,'.dz',Th,Vh,quadOrder);

cc = (uxec-uhxec).^2 + (uyec-uhyec).^2  + (uzec-uhzec).^2;
ErrH1 = integral3d(Th,cc,Vh,quadOrder);
ErrH1 = sqrt(ErrH1);