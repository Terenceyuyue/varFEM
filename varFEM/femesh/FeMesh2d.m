function Th = FeMesh2d(node,elem,bdStr)
%%FeMesh2d  mesh-related info in 2-D

if nargin==2, bdStr = [];  end

N = size(node,1); NT = size(elem,1);

%% Determine the number of boundary types
if ~iscell(bdStr) && ~isempty(bdStr)
    bdStr = {bdStr};
end
nbdType = length(bdStr)+1; % +1 for remaining boundary edges

%% edge
allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
totalEdge = sort(allEdge,2);
[i,j,s] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
edge = [j,i];  
NE = size(edge,1);

%% bdEdge (counterclockwise)
[~, i1, totalJ] = unique(totalEdge,'rows');
bdEdge = allEdge(i1(s==1),:);
bdEdgeIdx = find(s==1);
nbE = size(bdEdge,1);

%% Restore the counterclockwise orientation of edges on domain boundary
edge(bdEdgeIdx,:) = bdEdge;

%% Set up boundary edges
midbdEdge = (node(bdEdge(:,1),:) + node(bdEdge(:,2),:))/2;
x = midbdEdge(:,1); y = midbdEdge(:,2); %#ok<NASGU> 
bdEdgeType = cell(nbdType,1);  
bdEdgeIdxType = cell(nbdType,1);
bdNodeIdxType = cell(nbdType,1);
IdxRes = true(nbE,1); 
for i = 1:nbdType-1
    Idx = false(nbE,1);  % initialized as false at each iteration
    bdStri = bdStr{i};
    if mycontains(bdStri,'==')    % transform 'x==1' to 'abs(x-1)<1e-4'
        bd = bdStri;
        bd = strrep(bd,'==','-');
        if mycontains(bd,'|')
            bd = strrep(bd, '|', ')<1e-4 | abs(');
        end
        bdStri = ['abs(',   bd,   ')<1e-4'];
    end
    id = eval(bdStri);
    Idx(id) = true;   IdxRes(id) = false;
    bdEdgeType{i} = bdEdge(Idx,:);
    bdEdgeIdxType{i} = bdEdgeIdx(Idx);
    bdNodeIdxType{i} = unique(bdEdgeType{i});
end
bdEdgeType{end} = bdEdge(IdxRes,:);
bdEdgeIdxType{end} = bdEdgeIdx(IdxRes);
bdNodeIdxType{end} = unique(bdEdgeType{end});

%% elem2edge
elem2edge = reshape(totalJ,NT,3);

%% edge2elem
totalJelem = repmat((1:NT)',3,1);
i2(totalJ) = 1:length(totalJ); i2 = i2(:); % second occurence by overlapping 
edge2elem = totalJelem([i1,i2]);

%% ne
v12 = node(edge(:,1),:)-node(edge(:,2),:);
Ne = [-v12(:,2),v12(:,1)];
he = vecnorm(v12,2,2);
ne = Ne./he;

%% area
area = simplexvolume(node,elem);

%% Stored as a struct
Th.node = node; Th.elem = elem;
% boundary
Th.bdNodeIdx = unique(bdEdge);
Th.bdEdge = bdEdge; Th.bdEdgeIdx = bdEdgeIdx;
Th.bdEdgeType = bdEdgeType;  
Th.bdEdgeIdxType = bdEdgeIdxType;
Th.bdNodeIdxType = bdNodeIdxType;
% auxstructure
Th.edge = edge; 
Th.elem2edge = elem2edge;
Th.edge2elem = edge2elem;
% number
Th.N = N; Th.NT = NT; Th.NE = NE;
% geometric quantities
Th.ne = ne;  Th.he = he;  Th.Ne = Ne;
Th.diameter = max(he(elem2edge), [], 2);
Th.area = area;
% bdStr
Th.bdStr = bdStr;