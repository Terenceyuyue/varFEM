function w = Base2D(wStr,node,elem,Vh,quadOrder) %#ok<STOUT>
%% BASE2D returns base matrix for numerical integration
%
%  w = {w1,w2, ..., wn}, where n is the number of the basis functions
%  wi is a matrix of size NT*nG, where NT and nG are the numbers of
%  2D elements and quadrature points
%
% Copyright (C) Terence Yu.

if nargin == 3, Vh = 'P1'; quadOrder = 3; end  % default: P1
if nargin == 4, quadOrder = 3; end

%% Quadrature and gradbasis
wStr = lower(wStr); %#ok<NASGU> % lowercase string
NT = size(elem,1);

% Gauss quadrature rule
[lambda,weight] = quadpts(quadOrder);  %#ok<ASGLU>
nQuad = length(weight);  %#ok<NASGU>

% gradbasis
Dlambda = gradbasis(node,elem);
Dlambda1 = Dlambda(1:NT,:,1);
Dlambda2 = Dlambda(1:NT,:,2);
Dlambda3 = Dlambda(1:NT,:,3);
Dlambdax = [Dlambda1(:,1), Dlambda2(:,1), Dlambda3(:,1)]; %#ok<NASGU>
Dlambday = [Dlambda1(:,2), Dlambda2(:,2), Dlambda3(:,2)]; %#ok<NASGU>

%% P1-Lagrange
if  strcmpi(Vh, 'P1')
    Base2D_P1;
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    Base2D_P2;   % matlab script
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3')
    Base2D_P3; % matlab script
end

%% Morley
if  strcmpi(Vh, 'Morley')
    Base2D_Morley; % matlab script
end

%% Zienkiewicz
if  strcmpi(Vh, 'Zienkiewicz')
    Base2D_Zienkiewicz; % matlab script
end