function pde = NavierStokesdata(nu)

if nargin<1, nu = 1; end

% --------- given by the symbolic computation ------
[u1,u1x,u1y,u2,u2x,u2y,pe,px,py,f1,f2] = compute_rhs(nu);

% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = [u1(x,y), u2(x,y)];
    end
    function val = pexact(p)
        x = p(:,1); y = p(:,2);
        val = pe(x,y);
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = [f1(x,y), f2(x,y)];
    end
% Dirichlet boundary conditions
    function val = g_D(p)
        val = uexact(p);
    end
% Neumann boundary conditions ( right side hand )
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        val = [u1x(x,y), u1y(x,y),u2x(x,y), u2y(x,y)];
    end

    function val = Dpexact(p)
        x = p(:,1); y = p(:,2);
        val = [px(x,y), py(x,y)];
    end

pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du', @Du, ...
    'pexact', @pexact, 'Dpexact',@Dpexact, 'nu',nu);
end

function [u1,u1x,u1y,u2,u2x,u2y,p,px,py,f1,f2] = compute_rhs(nu)
syms x y;
% exact solution
phix = x.^2.*(x-1).^2;
phiy = y.^2.*(y-1).^2;
dphix = 2*x.*(x-1).*(2*x-1);
dphiy = 2*y.*(y-1).*(2*y-1);
u1 = phix.*dphiy;
u2 = -dphix.*phiy;
p = (2*x-1).*(2*y-1);
% derivative
u1x = diff(u1,x);      u1y = diff(u1,y);
u2x = diff(u2,x);      u2y = diff(u2,y);
px = diff(p,x);        py = diff(p,y);
u1xx = diff(u1x,x);    u1yy = diff(u1y,y);
u2xx = diff(u2x,x);    u2yy = diff(u2y,y);
% Lap u
Lapu1 = u1xx + u1yy;
Lapu2 = u2xx + u2yy;
Nu1 = u1*u1x + u2*u1y;
Nu2 = u1*u2x + u2*u2y;
% f
f1 = -nu*Lapu1 + Nu1 + px;
f2 = -nu*Lapu2 + Nu2 + py;
% f1 =  4*y - 4*nu*(2*y - 1)*(3*x^4 - 6*x^3 + 6*x^2*y^2 - 6*x^2*y + 3*x^2 ...
%     - 6*x*y^2 + 6*x*y + y^2 - y) + 8*x^3*y^2*(x - 1)^2*(2*x^2 - 3*x + 1)*(2*y^2 - 3*y + 1)^2 ...
%     - 4*x^3*y^2*(x - 1)^2*(y - 1)^2*(2*x^2 - 3*x + 1)*(6*y^2 - 6*y + 1) - 2;
% f2 = 4*x + 4*nu*(2*x - 1)*(6*x^2*y^2 - 6*x^2*y + x^2 - 6*x*y^2 + ...
%     6*x*y - x + 3*y^4 - 6*y^3 + 3*y^2) + 8*x^2*y^3*(y - 1)^2*(2*x^2 - 3*x + 1)^2*(2*y^2 - 3*y + 1)...
%     - 4*x^2*y^3*(x - 1)^2*(y - 1)^2*(6*x^2 - 6*x + 1)*(2*y^2 - 3*y + 1) - 2;


% convert to anonymous functions
u1 = matlabFunction(u1,'Vars',{x,y});
u2 = matlabFunction(u2,'Vars',{x,y});
p = matlabFunction(p,'Vars',{x,y});

u1x = matlabFunction(u1x,'Vars',{x,y}); u1y = matlabFunction(u1y,'Vars',{x,y});
u2x = matlabFunction(u2x,'Vars',{x,y}); u2y = matlabFunction(u2y,'Vars',{x,y});
px = matlabFunction(px,'Vars',{x,y}); py = matlabFunction(py,'Vars',{x,y});

f1 = matlabFunction(f1,'Vars',{x,y});
f2 = matlabFunction(f2,'Vars',{x,y});
end
