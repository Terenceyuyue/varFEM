function pde = Poissondata_afem()

% --------- given by the symbolic computation ------
[af,cf,u,ux,uy,rhs] = compute_rhs();

% coef1
    function val = a(p)
        x = p(:,1); y = p(:,2);
        val = af(x,y) + 0*x;
    end
% coef2
    function val = c(p)
        x = p(:,1); y = p(:,2);
        val = cf(x,y) + 0*x;
    end

% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = u(x,y);
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = rhs(x,y);
    end
% Dirichlet boundary conditions
    function val = g_D(p)
        val = uexact(p);
    end
% for Neumann boundary conditions
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        val = [ux(x,y), uy(x,y)];
    end

pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du',@Du, 'a',@a, 'c',@c);
end

function [a,c,u,ux,uy,f] = compute_rhs()
    syms x y;
    % coef
    a = 1 + 0*x;
    c = 0 + 0*x;
    % exact solution
    u = x*y*(1-x)*(1-y)*exp(-1000*((x-0.5)^2+(y-0.117)^2));
    %u = sin(pi*y)*sin(pi*x);
    % derivative
    ux = diff(u,x);  uy = diff(u,y);
    aux = a*ux;  auy = a*uy;
    % f
    f = -(diff(aux,x)+diff(auy,y)) + c*u;
    % convert to anonymous functions
    a = matlabFunction(a,'Vars',{x,y});
    c = matlabFunction(c,'Vars',{x,y}); 
    u = matlabFunction(u,'Vars',{x,y});    
    f = matlabFunction(f,'Vars',{x,y});    
    ux = matlabFunction(ux,'Vars',{x,y});  
    uy = matlabFunction(uy,'Vars',{x,y}); 
end