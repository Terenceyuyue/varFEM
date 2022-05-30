% Distmesh2d examples.
% Copyright (C) Terence YUE Yu.

clc;clear;close all

Domain = {@Circle_Domain,   @Rectangle_Domain,          @Ellipse_Domain,...
    @Polygon_Domain,  @Rectangle_Circle_Domain,   @Circle_Circle_Domain,...
    @Horn_Domain,     @Suspension_Domain,         @Upper_Circle_Domain, ...
    @Upper_Circle_Circle_Domain, @Wrench_Domain, @Square_sizefunc,...
    @NACA0012_airfoil, @Pie_hole, @Superellipse, @Implicit,...
    };
% very slow: Superellipse, Implicit
h = [
    0.20,   0.20,   0.20, ...
    0.10,   0.05,   0.05, ...
    0.01,   1.00,   0.20, ...
    0.10,   0.05,   0.01,...
    0.01,   0.005,...
    0.08,   0.75
    ];

for id = 1:14  % very slow: Superellipse, Implicit (id = 15,16)
    
    func = Domain{id};
    str = func2str(func);
    [fd,fh,BdBox,pfix] = feval(func);
    
    fprintf('\n ------ %s ------ \n', str);
    
    h0 = h(id);
    
    if mycontains(str, 'Polygon')
        [node,elem] = distmesh2d(fd,fh,h0,BdBox,pfix,pfix);
    else
        [node,elem] = distmesh2d(fd,fh,h0,BdBox,pfix);
    end
    
    clf, showmesh(node,elem)
    str1 = strrep(str,'_','\_');
    title(sprintf('%s',str1))
    pause(0.1)
    
    
end


