function w = varPlateBending(Th,pde,Vh,quadOrder)
%varPlateBending solves plate bending problem using Morley or Zienkiewicz element.
% 
%       -D_{ij} M_{ij}(w) + cw = f in \Omega,
%       Dirichlet boundary condition:
%               w = g1, grad(w)n = g2    on \Gamma.
%
%   References
%   K. Feng and Z.C. Shi. "Mathematical Theory of Elastic Structures", Springer-Verlag, 1996.
%
% Copyright (C)  Terence Yu.

% Quadrature orders for int1d and int2d
if nargin==2, Vh = 'Morley'; quadOrder = 3; end % default: Morley
if nargin==3, quadOrder = 3; end

para = pde.para; D = para.D; nu = para.nu; c = para.c;

%% Assemble stiffness matrix
% a(v,u) = D*(v11*u11+v22*u22+nu*(v22*u11+v11*u22)+2*(1-nu)*v12*u12)+cvu
% Omega
Coef  = {D, D, D*nu, D*nu, D*2*(1-nu), c};
Test  = {'v.dxx', 'v.dyy', 'v.dyy', 'v.dxx', 'v.dxy', 'v.val'};
Trial = {'u.dxx', 'u.dyy', 'u.dxx', 'u.dyy', 'u.dxy', 'u.val'};
kk = int2d(Th,Coef,Test,Trial,Vh,quadOrder);

%% Assemble the right hand side
% Omega
Coef = pde.f;  Test = 'v.val';
ff = int2d(Th,Coef,Test,[],Vh,quadOrder);

%% Dirichlet boundary conditions
bdStruct = Th.bdStruct;
if strcmpi(Vh,'Morley')    
    bdNodeIdx = bdStruct.bdNodeIdx; bdEdgeIdx = bdStruct.bdEdgeIdx;
    bdEdgeD = bdStruct.bdEdgeD;
    g_D = pde.g_D;  Dw = pde.Dw;
    
    N = Th.N; NE = Th.NE; NNdof = N+NE;
    fixedNode = [bdNodeIdx; bdEdgeIdx+N];
    isBdNode = false(NNdof,1); isBdNode(fixedNode) = true;
    freeDof = ~isBdNode;
    
    node = Th.node;
    nodeD = node(bdNodeIdx,:); wD = g_D(nodeD);
    z1 = node(bdEdgeD(:,1),:); z2 = node(bdEdgeD(:,2),:); zc = (z1+z2)./2;
    e = z1-z2;  % e = z2-z1
    ne = [-e(:,2),e(:,1)];  % scaled
    ne = ne./repmat(sqrt(sum(ne.^2,2)),1,2);
    wnD = sum(Dw(zc).*ne,2);
    
    w = zeros(NNdof,1); w(fixedNode) = [wD;wnD];
    ff = ff - kk*w;
end

if strcmpi(Vh,'Zienkiewicz')
    bdNodeIdx = bdStruct.bdNodeIdx;
    
    N = Th.N;
    id = [bdNodeIdx; bdNodeIdx+N; bdNodeIdx+2*N];
    isBdDof = false(3*N,1); isBdDof(id) = true;
    bdDof = isBdDof; freeDof = (~isBdDof);
    
    g_D = pde.g_D;  Dw = pde.Dw;
    nodeD = Th.node(bdNodeIdx,:);  
    wD = g_D(nodeD);  Dw = Dw(nodeD); %wxD = Dw(:,1); wyD = Dw(:,2);
    
    w = zeros(3*N,1); w(bdDof) = [wD;Dw(:)];
    ff = ff - kk*w;
end

%% Solver
w(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
