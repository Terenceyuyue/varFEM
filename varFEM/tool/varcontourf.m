function [c,h] = varcontourf(varargin)
%function [c,h] = varcontourf(x,y,node,elem,uh,n)
%function [c,h] = varcontourf(node,elem,uh,n)

[a,b] = size(varargin{1});
switch min(a,b)
    case 1  % the first form
        if nargin==5
            [x,y,node,elem,uh] = deal(varargin{:});
        end
        if nargin==6
            [x,y,node,elem,uh,n] = deal(varargin{:});
        end
    case 2
        if nargin==3
            [node,elem,uh] = deal(varargin{:});
        end
        if nargin==4
            [node,elem,uh,n] = deal(varargin{:});
        end
        xmin = min(node(:,1));  xmax = max(node(:,1));
        ymin = min(node(:,2));  ymax = max(node(:,2));
        x = linspace(xmin,xmax,401);
        y = linspace(ymin,ymax,401);
end

[x,y] = meshgrid(x,y);
uxy = varInterpolant2d(x,y,node,elem,uh); 
%uxy = pdeInterpolant2d(x,y,node,elem,uh); % it is not efficient
[c,h] = contourf(x,y,uxy,n); axis equal
color = jet(10);  colormap(color);