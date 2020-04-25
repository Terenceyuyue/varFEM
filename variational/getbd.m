function [bddof,bdval] = getbd(Th,g_D,Vh)
% g_D: Dirichlet function handle

if nargin==2, Vh = 'P1'; end % default: P1

% --------- Dirichlet boundary conditions ---------------
bdStruct = Th.bdStruct;
eD = bdStruct.eD; elemD = bdStruct.elemD; bdIndexD = bdStruct.bdIndexD;

node = Th.node;  N = size(node,1);
%% P1-Lagrange
if  strcmpi(Vh, 'P1') 
    NNdof = N;
    id = eD;
    isBdNode = false(NNdof,1); isBdNode(id) = true;
    bddof = (isBdNode); 
    pD = node(eD,:);
    bdval = g_D(pD);
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')
    NE = Th.NE; NNdof = N + NE;
    id = [eD; bdIndexD + N];
    isBdNode = false(NNdof,1); isBdNode(id) = true;
    bddof = (isBdNode); 
    z1 = node(elemD(:,1),:); z2 = node(elemD(:,2),:);   zc = (z1+z2)/2;  
    pD = node(eD,:);  
    uD = g_D(pD); uDc = g_D(zc); 
    bdval = [uD; uDc];
end

%% P3-Lagrange
if  strcmpi(Vh, 'P3')
    NE = Th.NE; NT = Th.NT; NNdof = N + 2*NE + NT;
    id = [eD; bdIndexD+N; bdIndexD+N+NE];
    isBdNode = false(NNdof,1); isBdNode(id) = true;
    bddof = (isBdNode); 
    z1 = node(elemD(:,1),:); z2 = node(elemD(:,2),:);  
    za = z1+(z2-z1)/3;  zb = z1+2*(z2-z1)/3;
    pD = node(eD,:);  
    uD = g_D(pD); uDa = g_D(za); uDb = g_D(zb);
    bdval = [uD; uDa; uDb];
end