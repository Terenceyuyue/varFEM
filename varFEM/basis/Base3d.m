function w = Base3d(wStr,node,elem3,Vh,quadOrder)
%% Base3d returns basis matrix for numerical integration
%
%  w = {w1,w2, ..., wn}, where n is the number of basis functions
%  wi is a matrix of size NT*ng, where NT and ng are the numbers of
%  3D elements and quadrature points
%

if nargin == 3, Vh = 'P1'; quadOrder = 4; end % default: P1
if nargin == 4, quadOrder = 4; end

wStr = lower(wStr); % lowercase string

%% P1-Lagrange
if  strcmpi(Vh, 'P1')
    w = Base3d_P1(wStr,node,elem3,quadOrder);
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    w = Base3d_P2(wStr,node,elem3,quadOrder);   
end