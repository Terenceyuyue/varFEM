function [fh,Nh] = interpEdge(f,Th,Vh)

node = Th.node;  elem1d = Th.elem1d;
N = size(node,1);
ncol = length(f(node(1,:)));

za = node(elem1d(:,1),:); zb = node(elem1d(:,2),:);
e = za-zb;
he = sqrt(sum(e.^2,2));
Ne = [-e(:,2), e(:,1)];
ne = Ne./repmat(he,1,2);

%% P1-Lagrange
if  strcmpi(Vh, 'P1')
    NNdof = N;

    fh = zeros(NNdof,ncol);    
    fh(elem1d(:,1),:) = f(za);
    fh(elem1d(:,2),:) = f(zb);
    
    Nh = zeros(NNdof,2);  % normal vectors
    Nh(elem1d(:,1),:) = ne;   
    Nh(elem1d(:,2),:) = ne;
    return;
end

%% P2-Lagrange
if  strcmpi(Vh, 'P2')    
    NNdof = Th.N + Th.NE;
    elem2dofEdge = [elem1d, Th.elem1dIdx+N];
     
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

%% P3-Lagrange
if  strcmpi(Vh, 'P3') 
    N = Th.N; NE = Th.NE; NT = Th.NT;
    NNdof = N + 2*NE + NT;
    elem2dofEdge = [elem1d, Th.elem1dIdx+N, Th.elem1dIdx+N+NE];
     
    fh = zeros(NNdof,ncol); 
    z1 = (2*za+zb)/3;
    z2 = (za+2*zb)/3;
    fh(elem2dofEdge(:,1),:) = f(za);
    fh(elem2dofEdge(:,2),:) = f(zb);
    fh(elem2dofEdge(:,3),:) = f(z1);
    fh(elem2dofEdge(:,4),:) = f(z2);
    
    Nh = zeros(NNdof,2);  % normal vectors
    Nh(elem2dofEdge(:,1),:) = ne;   
    Nh(elem2dofEdge(:,2),:) = ne;
    Nh(elem2dofEdge(:,3),:) = ne;
    Nh(elem2dofEdge(:,4),:) = ne;
    return;
end