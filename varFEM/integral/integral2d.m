function [I,elemI] = integral2d(Th,fh,Vh,quadOrder,varargin)
%  int_\Omega f(x,y) dxdy, where f(x,y) is apparoximated by the FEM
%  function:
%    f(x,y) \approx fh(x,y) = f1*Phi1(x,y) + f2*Phi2(x,y) + ....
%  Then,
%    I = f1 * \int_\Omega(Phi1) + f2 * \int_\Omega(Phi2) + ....
%
%  fh: - a function handle
%      - an FEM function
%      - Coef matrix
%

%I = assem2d(Th,fh,[],[],Vh,quadOrder);

wStr = '.val';
if ~isempty(varargin)
    wStr = varargin{1};
end

%% Preparation
node = Th.node; elem = Th.elem;  
NT = Th.NT;
% Sparse assembly index
[elem2dof,~,NNdof] = dof2d(Th,Vh);
% Guassian quadrature
[lambda,weight] = quadpts(quadOrder);
weight = weight(:)'; % must be a row vector
ng = length(weight);
% area
z1 = node(elem(:,1),:);
z2 = node(elem(:,2),:);
z3 = node(elem(:,3),:);
area = simplexvolume(node,elem);

%% Coef matrix
cf = fh;
% case 1: cf = uh is a finite element function
if isnumeric(cf) && length(cf)==NNdof
    phi = Base2d(wStr,node,elem,Vh,quadOrder); % phi1,phi2,phi3
    CoefMat = zeros(NT,ng);
    for ib = 1:length(phi)
        ui = cf(elem2dof(:,ib));
        CoefMat = CoefMat + repmat(ui,1,ng).*phi{ib};   % ui*phi_i
    end
end

% case 2: coef matrix
if isnumeric(cf) && size(cf,1)==NT && size(cf,2)==ng
    CoefMat = cf;
end

% case 3: function
if isnumeric(cf) && length(cf)==1 % a constant
    cf = @(pz) cf+0*pz(:,1);
end
if ~isnumeric(cf)
    CoefMat = zeros(NT,ng);
    for p = 1:ng
        pz = lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
        CoefMat(:,p) = cf(pz);
    end
end

%% Integral
% I = \int_\Omega Coef dx
%   = assem2d(Th,Coef,[],[],Vh,quadOrder)
% Weight matrix
ww = repmat(weight,NT,1);
elemI = area.*sum(ww.*CoefMat,2);
I = sum(elemI);