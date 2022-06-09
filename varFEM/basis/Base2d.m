function w = Base2d(wStr,node,elem,Vh,quadOrder) 
%% Base2D returns basis matrix for numerical integration
%
%  w = {w1,w2, ..., wn}, where n is the number of basis functions
%  wi is a matrix of size NT*ng, where NT and ng are the numbers of
%  2D elements and quadrature points
%

if nargin == 3, Vh = 'P1'; quadOrder = 3; end  % default: P1
if nargin == 4, quadOrder = 3; end

wStr = lower(wStr); % lowercase string

%% P1-Lagrange
if  strcmpi(Vh, 'P1')
    w = Base2D_P1(wStr,node,elem,quadOrder);
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    w = Base2D_P2(wStr,node,elem,quadOrder);   
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3')
    w = Base2D_P3(wStr,node,elem,quadOrder); 
end
% 
% %% Crouzeix-Raviart linear element
% if  strcmpi(Vh, 'CR')
%     Base2D_CR; % matlab script
% end
% 
% %% Morley
% if  strcmpi(Vh, 'Morley')
%     Base2D_Morley; % matlab script
% end
% 
% %% Zienkiewicz
% if  strcmpi(Vh, 'Zienkiewicz')
%     Base2D_Zienkiewicz; % matlab script
% end