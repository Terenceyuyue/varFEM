function bdStruct= setboundary(node,elem,varargin)
% Set type of boundary edges and returns structure of boundary information.
% varargin: string for Neumann boundary
%
% Copyright (C) Terence Yu.

if size(elem,2) == 4 % rectangle element
    elem = mat2cell(elem,ones(size(elem,1),1), 4);
end

%% Find boundary edges
% totalEdge
if ~iscell(elem)  % Triangulation
    allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
else
    shiftfun = @(verts) [verts(2:end),verts(1)];  % or shiftfun = @(verts) circshift(verts,-1);
    T1 = cellfun(shiftfun, elem, 'UniformOutput', false);
    v0 = horzcat(elem{:})'; % the starting points of edges
    v1 = horzcat(T1{:})'; % the ending points of edges
    allEdge = [v0,v1];
end
totalEdge = sort(allEdge,2);
% counterclockwise bdEdge
[~, i1,totalJ] = unique(totalEdge,'rows');     % first occurrence
i2(totalJ)= 1:length(totalJ); i2 = i2(:);
bdEdge = allEdge(i1(i1==i2),:);

%% Set up boundary edges
nE = size(bdEdge,1);
% initialize as Dirichlet (true for Dirichlet, false for Neumann)
Idx = true(nE,1);
nodebdEdge = (node(bdEdge(:,1),:) + node(bdEdge(:,2),:))/2;
x = nodebdEdge(:,1); y = nodebdEdge(:,2); %#ok<NASGU>
nvar = length(varargin); % 1 * size(varargin,2)
% nvar==0: all boundary edges are Dirichlet
if nvar==0
    Idx = true(nE,1); % true for Dirichlet
end
if nvar==1
    bdNeumann = varargin{1};
    % case 1: all Dirichlet
    if isempty(bdNeumann)
        Idx = true(nE,1);  % true for Dirichlet
    % case 2: all Neumann
    elseif strcmpi(bdNeumann,'all')
        Idx = false(nE,1);
    % case 3:partial Neumann
    else
        if mycontains(bdNeumann,'==')    % transform 'x==1' to 'abs(x-1)<1e-4'
            bdStr = bdNeumann;
            bdStr = strrep(bdStr,'==','-');
            if mycontains(bdStr,'|')
                bdStr = strrep(bdStr, '|', ')<1e-4 | abs(');
            end
            bdNeumann = ['abs(',   bdStr,   ')<1e-4'];
        end
        id = eval(bdNeumann);
        Idx(id) = false;
    end
end
if nvar>1
    for i = 1:nvar
        bdNeumann = varargin{i};
        if mycontains(bdNeumann,'==')    % transform 'x==1' to 'abs(x-1)<1e-4'
            bdStr = bdNeumann;
            bdStr = strrep(bdStr,'==','-');
            if mycontains(bdStr,'|')
                bdStr = strrep(bdStr, '|', ')<1e-4 | abs(');
            end
            bdNeumann = ['abs(',   bdStr,   ')<1e-4'];
        end
        id = eval(bdNeumann);
        Idx(id) = false;
    end
end

%% bdStruct
bdStruct.bdEdge = bdEdge;   % all boundary edges
bdStruct.bdEdgeD = bdEdge(Idx,:); % Dirichlet boundary edges
bdStruct.bdEdgeN = bdEdge(~Idx,:); % Neumann boundary edges
bdStruct.bdNodeIdx = unique(bdEdge); % all boundary nodes
bdStruct.bdNodeIdxD = unique(bdEdge(Idx,:)); % Dirichlet boundary nodes
bdStruct.bdNodeIdxN = unique(bdEdge(~Idx,:)); % Neumann boundary nodes
bdEdgeIdx = find(i1==i2);      % indices of all boundary edges
bdEdgeIdxD = bdEdgeIdx(Idx);  % indices of Dirichelt boundary edges
bdEdgeIdxN = bdEdgeIdx(~Idx); % indices of Neumann boundary edges
bdStruct.bdEdgeIdx = bdEdgeIdx;
bdStruct.bdEdgeIdxD = bdEdgeIdxD;
bdStruct.bdEdgeIdxN = bdEdgeIdxN;

%% Find boundary elements with Neumann boundary edge as the first side
% find boundary elements with boundary edge as the first side
if ~iscell(elem) % triangulation
    NT = size(elem,1);
    bdElemIdx = mod(i1(i1==i2),NT);
    bdElemIdx(bdElemIdx==0) = NT;
    elem2edge = reshape(totalJ ,NT,3);
    bdElem2edge = elem2edge(bdElemIdx,:);
    %idx1 = ((bdElem2edgeN(:,1)-bdEdgeIdx)==0);
    idx2 = ((bdElem2edge(:,2)-bdEdgeIdx)==0);
    idx3 = ((bdElem2edge(:,3)-bdEdgeIdx)==0);
    bdElem2edge(idx2 ,:) = bdElem2edge(idx2 ,[2,3,1]);
    bdElem2edge(idx3 ,:) = bdElem2edge(idx3 ,[3,1,2]);    
    
    bdStruct.bdElem2edge = bdElem2edge;
    bdStruct.bdElem2edgeD = bdElem2edge(Idx,:);
    bdStruct.bdElem2edgeN = bdElem2edge(~Idx,:);
end
