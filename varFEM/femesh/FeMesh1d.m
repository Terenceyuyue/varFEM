function Th = FeMesh1d(node,elem1d,bdStr)
%%FeMesh1d mesh-related info in 1-D
%

if nargin==2, bdStr = [];  end

%% Determine the number of boundary types
if ~iscell(bdStr) && ~isempty(bdStr)
    bdStr = {bdStr};
end
nbdType = length(bdStr)+1; % +1 for remaining boundary edges

%% end points
bdPoints = find(accumarray(elem1d(:),1) == 1);
nbP = size(bdPoints,1); 

%% Set up boundary points
nodePoints = node(bdPoints,:);
if size(node,2)==1 % 1-d
    x = nodePoints; %#ok<NASGU> 
end
if size(node,2)==2 % 2-d boundary edges
    x = nodePoints(:,1); %#ok<NASGU> 
    y = nodePoints(:,2); %#ok<NASGU> 
end
bdNodeIdxType = cell(nbdType,1);
IdxRes = true(nbP,1); 

for i = 1:nbdType-1
    Idx = false(nbP,1);  % initialized as false at each iteration
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
    bdNodeIdxType{i} = bdPoints(Idx,:);
end
if sum(~IdxRes)==nbP
    bdNodeIdxType(end) = [];    
else
    bdNodeIdxType{end} = bdPoints(IdxRes,:);
end

%% Stored as a struct
Th.node = node; Th.elem1d = elem1d; 
% end points
Th.bdNodeIdx = unique(bdPoints);
Th.bdNodeIdxType = bdNodeIdxType;
% numbers
Th.N = size(node,1); Th.nel = size(elem1d,1);
% bdStr
Th.bdStr = bdStr;