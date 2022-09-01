function pde = elasticitydata3d(para)

% ------ Lame constants ------
if nargin == 0
    lambda = 1;   mu = 1;
else
    lambda = para.lambda; mu = para.mu;
end

para = struct('lambda',lambda,'mu',mu);
[u1,u2,u3,f1,f2,f3,u1x,u1y,u1z,u2x,u2y,u2z,u3x,u3y,u3z,...
    sig11,sig22,sig33,sig12,sig13,sig23] = compute_rhs(para);

% rhs
    function val = f(p)  
        x = p(:,1); y = p(:,2);  z = p(:,3);
        val = [f1(x,y,z)+0*x, f2(x,y,z)+0*x, f3(x,y,z)+0*x];
    end
% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);  z = p(:,3);
        val = [u1(x,y,z)+0*x, u2(x,y,z)+0*x, u3(x,y,z)+0*x];
    end
% Dirichlet
    function val = g_D(p)
        val = uexact(p);
    end
% derivative of the exact solution
    function val =  Du(p)  % u = [u1,u2,u3]
        x = p(:,1); y = p(:,2);  z = p(:,3);
        val(:,1) = u1x(x,y,z)+0*x;  val(:,2) = u1y(x,y,z)+0*x;   val(:,3) = u1z(x,y,z)+0*x;
        val(:,4) = u2x(x,y,z)+0*x;  val(:,5) = u2y(x,y,z)+0*x;   val(:,6) = u2z(x,y,z)+0*x;
        val(:,7) = u3x(x,y,z)+0*x;  val(:,8) = u3y(x,y,z)+0*x;   val(:,9) = u3z(x,y,z)+0*x;
    end

    function val = g_N(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        val = [sig11(x,y,z)+0*x,sig22(x,y,z)+0*x,sig33(x,y,z)+0*x,...
            sig12(x,y,z)+0*x,sig13(x,y,z)+0*x,sig23(x,y,z)+0*x];
    end

pde = struct('lambda',lambda, 'mu',mu, 'f', @f, 'uexact',@uexact,...
    'g_D',@g_D,'Du',@Du,'g_N',@g_N);
end

function [u1,u2,u3,f1,f2,f3,u1x,u1y,u1z,u2x,u2y,u2z,u3x,u3y,u3z,...
    sig11,sig22,sig33,sig12,sig13,sig23] = compute_rhs(para)    
    lambda = para.lambda; mu = para.mu; 
    syms x y z;    
    
    % exact solution
    u1 = (-1+cos(2*pi*x))*sin(2*pi*y) + 1/(1+lambda)*sin(pi*x)*sin(pi*y);
    u2 = -(-1+cos(2*pi*y))*sin(2*pi*x) + 1/(1+lambda)*sin(pi*x)*sin(pi*y);
    u3 = -(-1+cos(2*pi*z))*sin(2*pi*x) + 1/(1+lambda)*sin(pi*x)*sin(pi*z);

    % derivative
    u1x = diff(u1,x);      u1y = diff(u1,y);      u1z = diff(u1,z);
    u2x = diff(u2,x);      u2y = diff(u2,y);      u2z = diff(u2,z);
    u3x = diff(u3,x);      u3y = diff(u3,y);      u3z = diff(u3,z);

    % epsilon_ij
    E11 = u1x;
    E22 = u2y;
    E33 = u3z;
    E12 = 0.5*(u1y + u2x);
    E13 = 0.5*(u1z + u3x);
    E23 = 0.5*(u2z + u3y);
    divu = E11 + E22 + E33;
    % sigma_ij
    sig11 = lambda*divu + 2*mu*E11;
    sig22 = lambda*divu + 2*mu*E22;
    sig33 = lambda*divu + 2*mu*E33;
    sig12 = 2*mu*E12;
    sig13 = 2*mu*E13;
    sig23 = 2*mu*E23;
    
    % f
    f1 = - (diff(sig11,x) + diff(sig12,y) + diff(sig13,z));
    f2 = - (diff(sig12,x) + diff(sig22,y) + diff(sig23,z));
    f3 = - (diff(sig13,x) + diff(sig23,y) + diff(sig33,z));
    
    % convert to anonymous functions
    u1 = matlabFunction(u1,'Vars',{x,y,z});    
    u2 = matlabFunction(u2,'Vars',{x,y,z});
    u3 = matlabFunction(u3,'Vars',{x,y,z});
    f1 = matlabFunction(f1,'Vars',{x,y,z});    
    f2 = matlabFunction(f2,'Vars',{x,y,z}); 
    f3 = matlabFunction(f3,'Vars',{x,y,z}); 
    u1x = matlabFunction(u1x,'Vars',{x,y,z});  
    u1y = matlabFunction(u1y,'Vars',{x,y,z});
    u1z = matlabFunction(u1z,'Vars',{x,y,z});
    u2x = matlabFunction(u2x,'Vars',{x,y,z});  
    u2y = matlabFunction(u2y,'Vars',{x,y,z}); 
    u2z = matlabFunction(u2z,'Vars',{x,y,z}); 
    u3x = matlabFunction(u3x,'Vars',{x,y,z});  
    u3y = matlabFunction(u3y,'Vars',{x,y,z}); 
    u3z = matlabFunction(u3z,'Vars',{x,y,z}); 
    sig11 = matlabFunction(sig11,'Vars',{x,y,z}); 
    sig22 = matlabFunction(sig22,'Vars',{x,y,z});
    sig33 = matlabFunction(sig33,'Vars',{x,y,z});
    sig12 = matlabFunction(sig12,'Vars',{x,y,z});
    sig13 = matlabFunction(sig13,'Vars',{x,y,z});
    sig23 = matlabFunction(sig23,'Vars',{x,y,z});
end