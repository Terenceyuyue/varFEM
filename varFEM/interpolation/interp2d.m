function fh = interp2d(f,Th,Vh)

node = Th.node;  elem = Th.elem;
N = size(node,1);
ncol = length(f(node(1,:)));

z1 = node(elem(:,1),:);
z2 = node(elem(:,2),:);
z3 = node(elem(:,3),:);

%% P1-Lagrange
if  strcmpi(Vh, 'P1')
    NNdof = N;

    fh = zeros(NNdof,ncol);    
    fh(elem(:), :) = [f(z1); f(z2); f(z3)];
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    % dof
    [elem2dof,~,NNdof] = dof2d(Th,Vh);
    % mid-points ( i-th vertice --- i-th edge )
    a1 = (z2+z3)/2;
    a2 = (z3+z1)/2;
    a3 = (z1+z2)/2;

    fh = zeros(NNdof,ncol);    
    fh(elem2dof(:), :) = [f(z1); f(z2); f(z3); f(a1); f(a2); f(a3)];
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3')
    % dof
    [elem2dof,~,NNdof] = dof2d(Th,Vh);
    % mid-points ( i-th vertice --- i-th edge )
    a1 = (2*z2+z3)/3;
    a2 = (2*z3+z1)/3;
    a3 = (2*z1+z2)/3;
    b1 = (z2+2*z3)/3;
    b2 = (z3+2*z1)/3;
    b3 = (z1+2*z2)/3;
    zc = (z1+z2+z3)/3;

    fh = zeros(NNdof,ncol);    
    fh(elem2dof(:), :) = [f(z1); f(z2); f(z3); ...
        f(a1); f(a2); f(a3); f(b1); f(b2); f(b3); f(zc)];
end