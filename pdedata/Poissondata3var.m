function pde = Poissondata3var
%
% Copyright (C) Terence Yu.

% --------- given by the symbolic computation ------
[af,cf,u,ux,uy,uz,rhs] = compute_rhs;

% coef1
    function val = a(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        val = af(x,y,z) + 0*x;
    end
% coef2
    function val = c(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        val = cf(x,y,z) + 0*x;
    end

% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        val = u(x,y,z);
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        val = rhs(x,y,z);
    end
% Dirichlet boundary conditions
    function val = g_D(p)
        val = uexact(p);
    end
% for Neumann boundary conditions
    function val = Du(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        val = [ux(x,y,z), uy(x,y,z), uz(x,y,z)];
    end

pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du',@Du, 'a',@a, 'c',@c);
end

function [a,c,u,ux,uy,uz,f] = compute_rhs()
    syms x y z;
    % coef
    a = 1 + 0*x; 
    c = 1 + 0*x;
    % exact solution
    u = cos(pi*x)*cos(pi*y)*cos(pi*z);  
    % derivative
    ux = diff(u,x);  uy = diff(u,y); uz = diff(u,z);
    aux = a*ux;  auy = a*uy; auz = a*uz;
    % f
    f = -(diff(aux,x)+diff(auy,y)+diff(auz,z)) + c*u;
    % convert to anonymous functions
    a = matlabFunction(a,'Vars',{x,y,z});
    c = matlabFunction(c,'Vars',{x,y,z}); 
    u = matlabFunction(u,'Vars',{x,y,z});    
    f = matlabFunction(f,'Vars',{x,y,z});    
    ux = matlabFunction(ux,'Vars',{x,y,z});  
    uy = matlabFunction(uy,'Vars',{x,y,z}); 
    uz = matlabFunction(uz,'Vars',{x,y,z}); 
end
