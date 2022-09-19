function w = Base2dBd(wStr,node,elem,Vh,quadOrder) 
%% Base2d returns basis matrix for numerical integration
%
%  w = {w1,w2, ..., wn}, where n is the number of basis functions
%  wi is a matrix of size NT*ng, where NT and ng are the numbers of
%  2D elements and quadrature points
%

if nargin == 3, Vh = 'P1'; quadOrder = 5; end  % default: P1
if nargin == 4, quadOrder = 5; end

wStr = lower(wStr); % lowercase string

%% P1-Lagrange
if  strcmpi(Vh, 'P1')
    w = Base2dBd_P1(wStr,node,elem,quadOrder);
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    w = Base2dBd_P2(wStr,node,elem,quadOrder);   
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3')
    w = Base2dBd_P3(wStr,node,elem,quadOrder); 
end