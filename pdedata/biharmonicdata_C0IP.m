function pde = biharmonicdata_C0IP()


% --------- given by the symbolic computation ------
[u,ux,uy,rhs] = compute_rhs();

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
% derivative
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        val = [ux(x,y), uy(x,y)];
    end


pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du', @Du);
end

function [u,ux,uy,f] = compute_rhs()
    syms x y;    
    u = 10*x^2*y^2*(1-x)^2*(1-y)^2*sin(pi*x);
    % derivative
    ux = diff(u,x);      uy = diff(u,y);
    Lapu = diff(u,x,2) + diff(u,y,2);
    LapLapu = diff(Lapu,x,2) + diff(Lapu,y,2);
    % f
    f = LapLapu;    
    % convert to anonymous functions
    u = matlabFunction(u,'Vars',{x,y});
    ux = matlabFunction(ux,'Vars',{x,y});
    uy = matlabFunction(uy,'Vars',{x,y});    
    f = matlabFunction(f,'Vars',{x,y});
end
