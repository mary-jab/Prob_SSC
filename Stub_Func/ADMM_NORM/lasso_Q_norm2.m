function [Z, history] = lasso_Q_nom2(A, Q, lambda0 , lambda1, rho, alpha, preZ)
% lasso  Solve lasso problem via ADMM
%
% [z, history] = lasso(A, b, lambda, rho, alpha);
%
% Solves the following problem via ADMM:
%
%   minimize 1/2*|| AC - b ||_2^2 + \lambda || C ||_1
%


t_start = tic;
% Global constants and defaults

QUIET    = 0;
MAC_ITER = 500;
ABSTOL   = 1e-6;
RELTOL   = 1e-4;

% Data preprocessing

[m, n] = size(A);
if nargin <7
    preZ = zeros(n,n);
end
% for id = 1:n
    
    b = A;
    % save a matriC-vector multiply
    Atb = A'*b;
    
    % ADMM solver
    
    C = zeros(n,n);
    z = preZ;%(:,id);%zeros(n,1);
    u = zeros(n,n);
    
    % cache the factorization  A == L*U
    [L, U] =  factor(A,Q, lambda1, rho);
    
    if ~QUIET
        fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
            'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
    end
    
    for k = 1:MAC_ITER
        % C-update
        q = Atb + rho*z - u;    % temporary value
        
%         if( m >= n )    % if skinny
            C = U \ (L \ q);
%         else            % if fat
%             C = q/rho - (A'*(U \ ( L \ (A*q) )))/rho^2;
%         end
        C(1:n+1:n*n)  = eps;
        
        % z-update with relaCation
        zold = z;
        C_hat = alpha*C + (1 - alpha)*zold;
        z = shrinkage(C_hat + u, (lambda0*ones(n,n)) /rho );% + lambda1*(1-Q)) /rho );
        
        % u-update
        u = u + (C_hat - z);
        
        % diagnostics, reporting, termination checks
        history.objval(k)  = objective(A, b, Q, lambda0, lambda1, C, z);%objective(A, b, lambda0, C, z);
        
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
    
    Z= z;
% end

if ~QUIET
    toc(t_start);
end

end

function p = objective(A, b, Q, lambda0, lambda1, C, z) %%%%%%%%%%%%%% need to be modified
p = 1/2*norm((A*C - b),2) + lambda0*norm(z,1);% +lambda1*norm((1-Q).*C,2) ;
end

function z = shrinkage(C, kappa)
z = max( 0, C - kappa ) - max( 0, -C - kappa );
end

% Cholesky factorization when A == L*U
function [L, U] = factor(A,Q, lambda2, rho)
[m, n] = size(A);
%% association matrix
aM =2*lambda2*(1-Q)'*(1-Q);

% if ( m >= n )    % if skinny
    L = chol( A'*A + aM +rho*speye(n), 'lower' );
% else            % if fat
%     L = chol( speye(m) + 1/rho*(A*A') , 'lower' );
% end

% force matlab to recognize the upper / lower triangular structure
L = sparse(L);
U = sparse(L');
end