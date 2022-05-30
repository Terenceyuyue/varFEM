function [fh,Nh] = interpEdge(f,Th,Vh)

node = Th.node;  elem1 = Th.elem1;
N = size(node,1);
ncol = length(f(node(1,:)));

za = node(elem1(:,1),:); zb = node(elem1(:,2),:);
e = za-zb;
he = sqrt(sum(e.^2,2));
Ne = [-e(:,2), e(:,1)];
ne = Ne./repmat(he,1,2);

%% P1-Lagrange
if  strcmpi(Vh, 'P1')
    NNdof = N;

    fh = zeros(NNdof,ncol);    
    fh(elem1(:,1),:) = f(za);
    fh(elem1(:,2),:) = f(zb);
    
    Nh = zeros(NNdof,2);  % normal vectors
    Nh(elem1(:,1),:) = ne;   
    Nh(elem1(:,2),:) = ne;
    return;
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')    
    NNdof = Th.N + Th.NE;
    elem2dofEdge = [elem1, Th.bdEdgeIdx1+N];
     
    fh = zeros(NNdof,ncol); 
    zc = (za+zb)/2;
    fh(elem2dofEdge(:,1),:) = f(za);
    fh(elem2dofEdge(:,2),:) = f(zb);
    fh(elem2dofEdge(:,3),:) = f(zc);
    
    Nh = zeros(NNdof,2);  % normal vectors
    Nh(elem2dofEdge(:,1),:) = ne;   
    Nh(elem2dofEdge(:,2),:) = ne;
    Nh(elem2dofEdge(:,3),:) = ne;
    return;
end

