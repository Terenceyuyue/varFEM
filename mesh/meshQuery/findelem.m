function findelem(node,elem,range)
%Findelem highlights some elements

hold on

% if size(elem,2) == 4
%     elem = mat2cell(elem,ones(size(elem,1),1), 4); % rectangle element
% end

if nargin==2
    range = (1:size(elem,1))';
end

center = zeros(length(range),2);
s = 1;
for iel = range(1):range(end)
    if iscell(elem)
        index = elem{iel};
    else
        index = elem(iel,:);
    end
    verts = node(index, :); verts1 = verts([2:end,1],:);
    area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
    area = 0.5*abs(sum(area_components));
    center(s,:) = sum((verts+verts1).*repmat(area_components,1,2))/(6*area);
    s = s+1;
end

marker = 'o';
if size(center,1)>40, marker = '.'; end
plot(center(:,1),center(:,2),marker,'LineWidth',1,'MarkerEdgeColor','b',...
    'MarkerFaceColor','y','MarkerSize',18);
text(center(:,1)-0.02,center(:,2),int2str(range),'FontSize',10,...
    'FontWeight','bold','Color','r');

hold off