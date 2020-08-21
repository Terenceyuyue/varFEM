function bdStruct = setboundary3(node,elem,varargin)
% setboundary3 sets type of boundary faces and returns structure of boundary
% information. 
% varargin: string for Neumann boundary
%
% Copyright (C) Terence Yu.

%% Find boundary faces
allFace = [elem(:,[2 4 3]);elem(:,[1 3 4]);elem(:,[1 4 2]);elem(:,[1 2 3])];
totalFace = sort(allFace,2);
[~,i1] = unique(totalFace,'rows');
[~,i2] = unique(totalFace(end:-1:1,:), 'rows');
i2 = size(totalFace,1)+1-i2;
bdFace = allFace(i1(i1==i2),:);  %counterclockwise bdFace

%% Set up boundary faces
nF = size(bdFace,1);
% initial as Dirichlet (true for Dirichlet, false for Neumann)
IdxD = true(nF,1);
nodebdFace = (node(bdFace(:,1),:) + node(bdFace(:,2),:) + node(bdFace(:,3),:))/3;
x = nodebdFace(:,1); y = nodebdFace(:,2); z = nodebdFace(:,3); %#ok<NASGU>
nvar = length(varargin); % 1 * size(varargin,2)
% note that length(varargin) = 1 for bdNeumann = [] or ''
if (nargin==2) || (~isempty(varargin{1})) 
    for i = 1:nvar 
        bdNeumann = varargin{i};
        id = eval(bdNeumann);
        IdxD(id) = false;
    end
end

bdStruct.bdFace = bdFace;   % all boundary faces
bdStruct.bdFaceD = bdFace(IdxD,:); % Dirichlet boundary faces
bdStruct.bdFaceN = bdFace(~IdxD,:); % Neumann boundary faces
bdFaceIdx = find(i1==i2);      % indices of all boundary faces
bdStruct.bdFaceIdx = bdFaceIdx;
bdStruct.bdFaceIdxD = bdFaceIdx(IdxD); % indices of Dirichelt boundary faces
bdStruct.bdFaceIdxN = bdFaceIdx(~IdxD); % indices of Neumann boundary faces

%% edge, elem2edge
% edge
allEdge = [elem(:,[1,2]); elem(:,[1,3]); elem(:,[1,4]);
           elem(:,[2,3]); elem(:,[2,4]); elem(:,[3,4])];
totalEdge = sort(allEdge,2);
[edge,~,totalJ] = unique(totalEdge,'rows');
NT = size(elem,1);
elem2edge = reshape(totalJ, NT, 6);
bdStruct.edge = edge;
bdStruct.elem2edge = elem2edge;

%% face2edge
allFace2edge = [elem2edge(:,[6,4,5]);
                elem2edge(:,[6,3,2]);
                elem2edge(:,[5,1,3]);
                elem2edge(:,[4,2,1])];
face2edge = allFace2edge(i1,:);
bdStruct.face2edge = face2edge;

%% Set up boundary edges
% initial as Dirichlet (true for Dirichlet, false for Neumann)
nE = size(edge,1);
Idx = true(nE,1);
nodeEdge = (node(edge(:,1),:) + node(edge(:,2),:))/2;
x = nodeEdge(:,1); y = nodeEdge(:,2); z = nodeEdge(:,3); %#ok<NASGU>
nvar = length(varargin); % 1 * size(varargin,2)
% note that length(varargin) = 1 for bdNeumann = [] or ''
if (nargin==2) || (~isempty(varargin{1})) 
    for i = 1:nvar 
        bdNeumann = varargin{i};
        id = eval(bdNeumann);
        Idx(id) = false;
    end
end

bdStruct.bdEdgeD = edge(Idx,:); % Dirichlet boundary edges
bdStruct.bdNodeDIdx = unique(edge(Idx,:)); % Dirichlet boundary nodes
bdStruct.bdEdgeDIdx = find(Idx);      % index of Dirichlet boundary edges
