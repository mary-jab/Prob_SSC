
function [Q, history] = findQ_lasso(Z, lambda, rho, alpha)
% lasso  Solve lasso problem via ADMM
%
% [z, history] = lasso(A, b, lambda, rho, alpha);
%
% Solves the following problem via ADMM:
%
%   minimize 1/2*|| Ax - b ||_2^2 + \lambda || x ||_1
%
% The solution is returned in the vector x.
%
% history is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% rho is the augmented Lagrangian parameter.
%
% alpha is the over-relaxation parameter (typical values for alpha are
% between 1.0 and 1.8).
%
%
% More information can be found in the paper linked at:
% http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%

t_start = tic;

QUIET    = 1;
MAX_ITER = 1000;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;
A = Z;

[m, n] = size(A);

% save a matrix-vector multiply

for it = 1:n
    zi = Z(:,it);
    I = ones(m,1);
    Atb = zi*zi'*I;
    
    qi = zeros(m,1);
    c = zeros(m,1);
    u = zeros(m,1);
    
    % cache the factorization
    [L U] = factor(zi, rho);
    
    if ~QUIET
        fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
            'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
    end
    
    for k = 1:MAX_ITER
        
        % qi-update
        q = Atb + rho*(c - u);    % temporary value
        if( m >= n )    % if skinny
            qi = U \ (L \ q);
        else            % if fat
            qi = q/rho - (zi*(U \ ( L \ (zi'*q) )))/rho^2;
        end
%         qi = pos(qi);
        qi = abs(qi);
        qi = qi/max(qi);
        % z-update with relaxation
        zold = c;
        qi_hat = alpha*qi + (1 - alpha)*zold;
        c = shrinkage(qi_hat + u, lambda/rho);
        
        % u-update
        u = u + (qi_hat - c);
        
        % diagnostics, reporting, termination checks
        history.objval(k)  = objective(zi, lambda, qi, c);
        
        history.r_norm(k)  = norm(qi - c);
        history.s_norm(k)  = norm(-rho*(c - zold));
        
        history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(qi), norm(-c));
        history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
        
        if ~QUIET
            fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
                history.r_norm(k), history.eps_pri(k), ...
                history.s_norm(k), history.eps_dual(k), history.objval(k));
        end
        
        if (history.r_norm(k) < history.eps_pri(k) && ...
                history.s_norm(k) < history.eps_dual(k))
            break;
        end
        
    end
    Q(:,it) = c;
end
if ~QUIET
    toc(t_start);
end
end

function p = objective(A, lambda, x, z)
p = ( 1/2*sum(((1-x)'*A).^2) + lambda*norm(z,1) );
end

function z = shrinkage(x, kappa)
z = max( 0, x - kappa ) - max( 0, -x - kappa );
end

function [L U] = factor(A, rho)
[m, n] = size(A);
if ( m >= n )    % if skinny
    L = chol( A'*A + rho*speye(n), 'lower' );
else            % if fat
    L = chol( speye(m) + 1/rho*(A*A'), 'lower' );
end

% force matlab to recognize the upper / lower triangular structure
L = sparse(L);
U = sparse(L');
end

