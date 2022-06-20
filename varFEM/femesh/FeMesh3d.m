function Th = FeMesh3d(node,elem3,bdStr)
%%FeMesh3d  mesh-related info in 3-D

if nargin==2, bdStr = []; end

% node, elem
[elem3fix,Idx] = fixorder3(node,elem3);  % with reverse normal vector
if length(Idx)>1  % only for uniformrefine3
    elem3 = elem3fix;
    elem3(:,[2 3]) = elem3(:,[3 2]);
end

% % WARNING: the output of fixorder3.m can not be used in uniformrefine3.m
Th.node = node; Th.elem3 = elem3;

%% Determine the number of boundary types
if ~iscell(bdStr) && ~isempty(bdStr)
    bdStr = {bdStr};
end
nbdType = length(bdStr)+1; % +1 for remaining boundary faces

%% bdFace (counterclockwise)
allFace = [elem3(:,[2 4 3]);elem3(:,[1 3 4]);elem3(:,[1 4 2]);elem3(:,[1 2 3])];
totalFace = sort(allFace,2);
[~,i1,totalJ] = unique(totalFace,'rows');
i2(totalJ) = 1:length(totalJ); i2 = i2(:); % second occurence by overlapping
bdFace = allFace(i1(i1==i2),:);
bdFaceIdx = find(i1==i2);      % index of all boundary faces

%% Set up boundary faces
nF = size(bdFace,1);
% initial as Dirichlet (true for Dirichlet, false for Neumann)
midbdFace = (node(bdFace(:,1),:) + node(bdFace(:,2),:) + node(bdFace(:,3),:))/3;
x = midbdFace(:,1); y = midbdFace(:,2); z = midbdFace(:,3); %#ok<NASGU>
bdFaceType = cell(nbdType,1);
bdFaceIdxType = cell(nbdType,1);
bdNodeIdxType = cell(nbdType,1);
IdxRes = true(nF,1);
for i = 1:nbdType-1
    Idx = false(nF,1);  % initialized as false at each iteration
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
    bdFaceType{i} = bdFace(Idx,:);
    bdFaceIdxType{i} = bdFaceIdx(Idx);
    bdNodeIdxType{i} = unique(bdFaceType{i});
end
bdFaceType{end} = bdFace(IdxRes,:);
bdFaceIdxType{end} = bdFaceIdx(IdxRes);
bdNodeIdxType{end} = unique(bdFaceType{end});

%% edge, elem2edge
% edge
allEdge = [elem3(:,[1,2]); elem3(:,[1,3]); elem3(:,[1,4]);
    elem3(:,[2,3]); elem3(:,[2,4]); elem3(:,[3,4])];
totalEdge = sort(allEdge,2);
[edge,~,totalJE] = unique(totalEdge,'rows');
NT = size(elem3,1);
elem2edge = reshape(totalJE, NT, 6);

%% face2edge
allFace2edge = [elem2edge(:,[6,4,5]);  % face1
    elem2edge(:,[6,3,2]);  % face2
    elem2edge(:,[5,1,3]);  % face3
    elem2edge(:,[4,2,1])]; % face4
face2edge = allFace2edge(i1,:);
bdFace2edgeType = cell(nbdType,1);
for s = 1:nbdType
    bdFace2edgeType{s} = face2edge(bdFaceIdxType{s},:);
end


%% Stored as a struct
% auxstructure
Th.edge = edge;
Th.elem2edge = elem2edge;
Th.face2edge = face2edge;
% boundary face
Th.bdFaceType = bdFaceType;
Th.bdFaceIdxType = bdFaceIdxType;
Th.bdFace2edgeType = bdFace2edgeType;
Th.bdNodeIdxType = bdNodeIdxType;
% boundary edge
bdEdgeIdxType = cell(nbdType,1);
bdEdgeType = cell(nbdType,1);
for s = 1:nbdType
    bdEdgeIdxType{s} = unique(bdFace2edgeType{s});
    bdEdgeType{s} = edge(bdEdgeIdxType{s},:);
end
Th.bdEdgeIdxType = bdEdgeIdxType;
Th.bdEdgeType = bdEdgeType;
% numbers
Th.N = size(node,1); Th.NT = size(elem3,1);  Th.NE = size(edge,1);
% bdStr
Th.bdStr = bdStr;