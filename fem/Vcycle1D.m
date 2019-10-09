% ---------------------------------------------------------------------
%                         Vcycle1D.m 
% ---------------------------------------------------------------------
function e = Vcycle1D(A,r,Ai,Pro,Res)
% Solve the residual equation Ae = r by multigrid V-cycle method

J = length(Ai); % level length

% If the problem is small enough, solve it directly
if J<=1
    e = A\r;    return;
end

ri = cell(J,1);            % residual in each level
ei = cell(J,1);            % correction in each level
ri{J} = r;

% ---------- Correction in each level ----------
option = 'forward';
for j = J:-1:2
    % % pre-smoothing: one step
    % ei{j} = R{j}*ri{j};
    ei{j} = smoother(Ai{j},ri{j},option);
    % update and restrict residual
    ri{j-1} = Res{j-1}*(ri{j}-Ai{j}*ei{j});
end
ei{1} = Ai{1}\ri{1}; % exact solver in the coarsest level

% -------------- prolongation and correction ------------
option = 'backward';
for j = 2:J
    % prolongation and correction
    ei{j} = ei{j}+Pro{j-1}*ei{j-1};
    % % post-smoothing: one step
    % ei{j} = ei{j} + R{j}'*(ri{j}-Ai{j}*ei{j});
    rij = ri{j}-Ai{j}*ei{j};
    ei{j} = ei{j} + smoother(Ai{j},rij,option);
end
e = ei{J};
end % end of Vcycle1D