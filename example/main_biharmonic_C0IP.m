%%% biharmonic_C0IP: Implementation of the symmetric quadratic C0 interior 
% penalty method for the biharmonic equation: \Delta^2 u = f
% Ref: Bringmann, P., Carstensen, C. and Streitberger, J., Parameter-free 
% implementation of the quadratic C0 interior penalty method for the biharmonic
% equation, arXiv:2209.05221v1

clc;clear;close all
tic;
%% Mesh
[node,elem] = squaremesh([0 1 0 1],0.1);
bdStr = [];
Th = FeMesh2d(node,elem,bdStr);

%% PDE data
pde = biharmonicdata_C0IP();

%% Finite element space
Vh = 'P2'; quadOrder = 5;

%% Piecewise biliner form
Coef  = {1,1,1,1};
Test  = {'v.dxx','v.dxy','v.dyx','v.dyy'};
Trial = {'u.dxx','u.dxy','u.dyx','u.dyy'};
A = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);

%% Preparations for jump and penalty terms
% numbers
N = Th.N; NT = Th.NT; NE = Th.NE;  NNdof = N + NE;  Ndofe = 9;
% Gauss quadrature in 1-D 
[~,w1d] = quadpts1(quadOrder);  ng = length(w1d);
ww = repmat(w1d,NE,1);
% auxiliary 
bdEdgeIdx = Th.bdEdgeIdx;
sgnedge = ones(NE,1);
sgnedge(bdEdgeIdx) = sign(Th.bdEdge(:,2)-Th.bdEdge(:,1));
he = Th.he;  
ne = repmat(sgnedge,1,2).*Th.ne;
nx = ne(:,1);  ny = ne(:,2);
nxg = repmat(nx,1,ng);  nyg = repmat(ny,1,ng);
% some parameters
a = 4;  % for stability
ce = 0.5*ones(NE,1);  ce(bdEdgeIdx) = 1;
se = 3*a/4*2*ce.*he.^2.*sum(1./Th.area(Th.edge2elem(:,1:2)),2); % penalty
Cav = ones(NE,1);  Cav(bdEdgeIdx) = 2;  % adjustment for average
Cavg = repmat(Cav,1,ng);  

%% edgeBase: left and right evaluations of basis functions on macro elements
wStr = {'.dx', '.dy', '.dxx', '.dxy', '.dyy'};
[edgeBaseM,edgeBaseP,edge2dof] = edgeBase(wStr,Th,Vh,quadOrder);
[wxM,wyM,wxxM,wxyM,wyyM] = deal(edgeBaseM{:});
[wxP,wyP,wxxP,wxyP,wyyP] = deal(edgeBaseP{:});

%% Sparse assembly index
ii = reshape(repmat(edge2dof, Ndofe,1), [], 1);
jj = repmat(edge2dof(:), Ndofe, 1);

%% Stiffness matrices J and C for jump and penalty terms
KJ = zeros(NE,Ndofe^2);  KC = zeros(NE,Ndofe^2);
s = 1;
for i = 1:Ndofe
    for j = 1:Ndofe
        % J
        % J: first term 
        v1i = wxM{i}-wxP{i};   
        u1j = 0.5*Cavg.*(wxxM{j}+wxxP{j}).*nxg ...
            + 0.5*Cavg.*(wxyM{j}+wxyP{j}).*nyg;
        % J: second term
        v2i = wyM{i}-wyP{i};   
        u2j = 0.5*Cavg.*(wxyM{j}+wxyP{j}).*nxg ...
            + 0.5*Cavg.*(wyyM{j}+wyyP{j}).*nyg;
        % J: combined
        KJ(:,s) = he.*sum(ww.*(v1i.*u1j + v2i.*u2j), 2);

        % C
        vi = (wxM{i}-wxP{i}).*nxg + (wyM{i}-wyP{i}).*nyg;
        uj = (wxM{j}-wxP{j}).*nxg + (wyM{j}-wyP{j}).*nyg;
        KC(:,s) = se.*sum(ww.*vi.*uj, 2);

        s = s+1;
    end
end

%% Assemble the jump and penalty terms
J = sparse(ii,jj,KJ(:),NNdof,NNdof);
C = sparse(ii,jj,KC(:),NNdof,NNdof);
kk = A - J - J' + C;

%% Assemble the right-hand side
Coef = pde.f;  Test = 'v.val';
ff = assem2d(Th,Coef,Test,[],Vh,quadOrder);

%% Apply Dirichlet boundary conditions
g_D = pde.g_D;
on = 1;
uh = apply2d(on,Th,kk,ff,Vh,g_D);

%% Plot
figure,
showresult(node,elem,pde.uexact,uh);
toc