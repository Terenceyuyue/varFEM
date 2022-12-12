function eta = Poisson_indicator(Th,uh,pde,Vh,quadOrder)
%%This function returns the local error indicator of Poisson equation with
% homogeneous Dirichlet boundary condition in 2-D.
%

%% preparation for the computation
elem2edge = Th.elem2edge;
he = Th.he;  
diameter = Th.diameter;
NT = Th.NT;

%% elementwise residuals
fc = interp2dMat(pde.f,'.val',Th,Vh,quadOrder);
uxxc = interp2dMat(uh,'.dxx',Th,Vh,quadOrder);
uyyc = interp2dMat(uh,'.dyy',Th,Vh,quadOrder);
Coef = (fc + uxxc + uyyc).^2;
[~,elemIh] = integral2d(Th,Coef,Vh,quadOrder);
elemRes = diameter.^2.*elemIh;

%% elementwise interior and exterior evaluations at quadrature points
[elemuhxM,elemuhxP,elemnx,elemny] = elem2edgeInterp('.dx',Th,uh,Vh,quadOrder);
[elemuhyM,elemuhyP] = elem2edgeInterp('.dy',Th,uh,Vh,quadOrder);

%% elementwise evaluations of the jump integral
elem2Jumpx = elemuhxM - elemuhxP;
elem2Jumpy = elemuhyM - elemuhyP;
[~,weight1d] = quadpts1(quadOrder);
ng = length(weight1d);
% elemJump = zeros(NT,1);
% for i = 1:3 % loop of triangle sides
%     hei = he(elem2edge(:,i));
%     id = (1:ng)+(i-1)*ng;
%     cei = hei;
%     neix = elemnx(:,id); 
%     neiy = elemny(:,id);
%     Jumpnx = elem2Jumpx(:,id).*neix;
%     Jumpny = elem2Jumpy(:,id).*neiy;
%     Jumpn = (Jumpnx+Jumpny).^2;    
%     elemJump = elemJump + cei.*hei.*(Jumpn*weight1d(:));
% end
id3 = elem2edge(:, kron([1 2 3], ones(1,ng)));
he3 = he(id3);
ce3 = he3;
ww3 = repmat(weight1d,NT,3);
Jumpnx = elem2Jumpx.*elemnx;
Jumpny = elem2Jumpy.*elemny;
Jumpn = (Jumpnx+Jumpny).^2;
elemJump = sum(ww3.*ce3.*he3.*Jumpn,2);

%% Local error indicator
eta = (abs(elemRes) + elemJump).^(1/2);