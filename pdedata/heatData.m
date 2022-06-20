function pde = heatData

% --------- given by the symbolic computation ------
[u,ux,uy,rhs] = compute_rhs;

% exact solution
    function val = uexact(p,t0)
        x = p(:,1); y = p(:,2); t = t0 + 0*x;
        val = u(x,y,t);
    end
% right side hand function
    function val = f(p,t0)
        x = p(:,1); y = p(:,2); t = t0 + 0*x;
        val = rhs(x,y,t);
    end
% Dirichlet boundary conditions
    function val = g_D(p,t0)
        val = uexact(p,t0);
    end
% for Neumann boundary conditions
    function val = Du(p,t0)
        x = p(:,1); y = p(:,2); t = t0 + 0*x;
        val = [ux(x,y,t), uy(x,y,t)];
    end

pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du',@Du);
end

function [u,ux,uy,f] = compute_rhs()
    syms x y t;
  
    % exact solution
    %u = exp(-x-y-t);
    u = sin(pi*x)*sin(y)*exp(-t);     
    % derivative
    ut = diff(u,t);
    ux = diff(u,x);  uy = diff(u,y);
    Lapu = diff(ux,x) + diff(uy,y);
    % f
    f = ut-Lapu;

    % convert to anonymous functions
    u = matlabFunction(u,'Vars',{x,y,t});    
    f = matlabFunction(f,'Vars',{x,y,t});    
    ux = matlabFunction(ux,'Vars',{x,y,t});  
    uy = matlabFunction(uy,'Vars',{x,y,t}); 
end
