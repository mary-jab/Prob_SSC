function [z, history] = lasso(A, b, lambda, rho, alpha, id)
% lasso  Solve lasso problem via ADMM
%
% [z, history] = lasso(A, b, lambda, rho, alpha);
%
% Solves the following problem via ADMM:
%
%   minimize 1/2*|| AC - b ||_2^2 + \lambda || C ||_1
%
% The solution is returned in the vector C.
%
% history is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% rho is the augmented Lagrangian parameter.
%
% alpha is the over-relaCation parameter (typical values for alpha are
% between 1.0 and 1.8).
%
%
% More information can be found in the paper linked at:
% http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%

t_start = tic;
% Global constants and defaults

QUIET    = 0;
MAC_ITER = 1000;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;

% Data preprocessing

[m, n] = size(A);

% save a matriC-vector multiply
Atb = A'*b;
% ADMM solver

C = zeros(n,1);
z = zeros(n,1);
u = zeros(n,1);

% cache the factorization
[L U] = factor(A, rho);

if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

for k = 1:MAC_ITER

    % C-update
    q = Atb + rho*(z - u);    % temporary value
    if( m >= n )    % if skinny
       C = U \ (L \ q);
    else            % if fat
       C = q/rho - (A'*(U \ ( L \ (A*q) )))/rho^2;
    end
     C(id) = eps;
%     C = C - diag(diag(C));

    % z-update with relaCation
    zold = z;
    C_hat = alpha*C + (1 - alpha)*zold;
    z = shrinkage(C_hat + u, lambda/rho);

    % u-update
    u = u + (C_hat - z);

    % diagnostics, reporting, termination checks
    history.objval(k)  = objective(A, b, lambda, C, z);

    history.r_norm(k)  = norm(C - z);
    history.s_norm(k)  = norm(-rho*(z - zold));

    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(C), norm(-z));
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

if ~QUIET
    toc(t_start);
end
end

function p = objective(A, b, lambda, C, z)
    p = ( 1/2*sum((A*C - b).^2) + lambda*norm(z,1) );
end

function z = shrinkage(C, kappa)
    z = max( 0, C - kappa ) - max( 0, -C - kappa );
end

function [L, U] = factor(A, rho)
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