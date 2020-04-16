function u = Applyboundary2Dvec(Th,kk,ff,pde,stru,feSpace)

if nargin==4, stru = 'u.val'; feSpace = 'P1'; end
if nargin==5, feSpace = 'P1'; end % default: P1

node = Th.node;  N = size(node,1);

% --------- Dirichlet boundary conditions ---------------
bdStruct = Th.bdStruct;
eD = bdStruct.eD; elemD = bdStruct.elemD; bdIndexD = bdStruct.bdIndexD;
g_D = pde.g_D;

NNdof = size(ff,1)/2;

%% P1-Lagrange
if  strcmpi(feSpace, 'P1')
    id0 = eD;
    if strcmpi(stru, 'u1.val'), id = id0; end
    if strcmpi(stru, 'u2.val'), id = id0; id = id + NNdof; end    
    if strcmpi(stru, 'u.val'),  id = id0; id = [id; id + NNdof]; end
    isBdNode = false(2*NNdof,1); isBdNode(id) = true;
    bdDof = (isBdNode); freeDof = (~isBdNode);
    pD = node(eD,:);
    u = zeros(2*NNdof,1); uD = g_D(pD); u(bdDof) = uD(:);
end

%% P2-Lagrange
if  strcmpi(feSpace, 'P2')
    id0 = [eD; bdIndexD + N];
    if strcmpi(stru, 'u1.val'), id = id0; end
    if strcmpi(stru, 'u2.val'), id = id0; id = id + NNdof; end    
    if strcmpi(stru, 'u.val'),  id = id0; id = [id; id + NNdof]; end
    isBdNode = false(2*NNdof,1); isBdNode(id) = true;
    bdDof = (isBdNode); freeDof = (~isBdNode);
    z1 = node(elemD(:,1),:); z2 = node(elemD(:,2),:);   zc = (z1+z2)/2;  
    pD = node(eD,:);  
    uD = g_D(pD); uDc = g_D(zc);  uDvec = [uD; uDc];
    u = zeros(2*NNdof,1); u(bdDof) = uDvec(:);
end

%% P3-Lagrange
if  strcmpi(feSpace, 'P3')
    NE = Th.NE;
    id0 = [eD; bdIndexD+N; bdIndexD+N+NE];
    if strcmpi(stru, 'u1.val'), id = id0; end
    if strcmpi(stru, 'u2.val'), id = id0; id = id + NNdof; end    
    if strcmpi(stru, 'u.val'),  id = id0; id = [id; id + NNdof]; end
    isBdNode = false(2*NNdof,1); isBdNode(id) = true;
    bdDof = (isBdNode); freeDof = (~isBdNode);
    z1 = node(elemD(:,1),:); z2 = node(elemD(:,2),:);  
    za = z1+(z2-z1)/3;  zb = z1+2*(z2-z1)/3;
    pD = node(eD,:);  
    uD = g_D(pD); uDa = g_D(za); uDb = g_D(zb);
    uDvec = [uD; uDa; uDb]; 
    u = zeros(2*NNdof,1); u(bdDof) = uDvec(:);
end

ff = ff - kk*u;

%% ------------------ Solver ----------------
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
