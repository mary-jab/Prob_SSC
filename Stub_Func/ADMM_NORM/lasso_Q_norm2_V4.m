function [Z, history] = lasso_Q_nom2(A, Q, lambda0 , lambda1, rho, alpha, preZ, options)
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

QUIET    = 1;
MAC_ITER = 200;
ABSTOL   = 1e-6;
RELTOL   = 1e-4;

% Data preprocessing

[D, N] = size(A);
if nargin <7
    preZ = zeros(N,N);
end
% for id = 1:n

b = A;
% save a matriC-vector multiply

% ADMM solver
if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
        'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

thr1 = 2*10^-4;
thr2 = 2*10^-4;

if (~options.ouliers)
    Atb = A'*b;
    
    C = zeros(N,N);
    z = preZ;%(:,id);%zeros(n,1);
    u = zeros(N,N);
    oulier = 0;
    % cache the factorization  A == L*U
    [L, U] =  factor(A,Q, lambda1, rho);
    for k = 1:MAC_ITER
        % C-update
        q = Atb + rho*z - u;    % temporary value
        
        %         if( m >= n )    % if skinny
        C = U \ (L \ q);
        %         else            % if fat
        %             C = q/rho - (A'*(U \ ( L \ (A*q) )))/rho^2;
        %         end
        C = C - diag(diag(C));
        C_hat = C;
        
        if (options.affine)
            C_hat = alpha*C + (1 - alpha)*z;
        end
        % z-update with relaCation
        zold = z;
        z = shrinkage(C_hat + u, (lambda0*ones(N,N)) /rho );% + lambda1*(1-Q)) /rho );
        z = z - diag(diag(z));
        
        % u-update
        u = u + (C_hat - z);
        
        % diagnostics, reporting, termination checks
        err1(k) = errorCoef(z,C);
        err2(k) = errorLinSys(A,z);
        
        if ~QUIET
            fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n', k, ...
                err1(k), thr1, ...
                err2(k), thr2);
        end
        
        if ((err1(k) < thr1 && err2(k) < thr2))
            break;
        end
        
    end
else
    gamma = alpha / norm(A,1);
    P = [A eye(D)/gamma];
    Qext = [Q eye(N,D)/gamma ];
    Atb = P'*b;
    
    C = zeros(N,N);
    z = [preZ; zeros(D, N)];%(:,id);%zeros(n,1);
    u = zeros(N+D,N);
    oulier = 0;
    % cache the factorization  A == L*U
    [L, U] =  factor(P,Qext, lambda1, rho);
    
    for k = 1:MAC_ITER
        % C-update
        q = Atb + rho*z - u;
        C = U \ (L \ q);
        C(1:N,:)  = C(1:N,:) - diag(diag(C(1:N,:)));
        C_hat = C;
        if (options.affine)
            C_hat = alpha*C + (1 - alpha)*z;
        end
        % z-update with relaCation
        zold = z;
        z = shrinkage(C_hat + u, (lambda0*ones(N+D,N)) /rho );% + lambda1*(1-Q)) /rho );
        z(1:N,:) = z(1:N,:) - diag(diag(z(1:N,:)));
        % u-update
        u = u + (C_hat - z);
        
        history.objval(k)  = objective(P, b, Qext, lambda0, lambda1, C, z);%objective(A, b, lambda0, C, z);
        history.r_norm(k)  = norm(C - z);
        history.s_norm(k)  = norm(-rho*(z - zold));
        
        history.eps_pri(k) = sqrt(N)*ABSTOL + RELTOL*max(norm(C), norm(-z));
        history.eps_dual(k)= sqrt(N)*ABSTOL + RELTOL*norm(rho*u);
        
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
    
end


Z= z(1:N,:);
% end
    fprintf('err1: %2.4f, err2: %2.4f, iter: %3.0f \n',err1(end),err2(end),k);

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

%  if ( m >= n )    % if skinny
L = chol( A'*A + aM +rho*speye(n), 'lower' );
%  else            % if fat
%      L = chol( speye(m) + 1/rho*((A*A') +aM) , 'lower' );
%  end

% force matlab to recognize the upper / lower triangular structure
L = sparse(L);
U = sparse(L');
end


function err = errorCoef(Z,C)
err = max(max( abs(Z-C) ));
end

function err = errorLinSys(P,Z)
[R,N] = size(Z);
if (R > N)
    E = P(:,N+1:end) * Z(N+1:end,:);
    Y = P(:,1:N);
    Y0 = Y - E;
    C = Z(1:N,:);
else
    Y = P;
    Y0 = P;
    C = Z;
end

for i = 1:size(Y0,2)
    n(i) = norm(Y0(:,i));
    Yn(:,i) = Y0(:,i) ./ n(i);
end

M = repmat(n,size(Y,1),1);
S = Yn - Y * C ./ M;
err = sqrt( max( sum( S.^2,1 ) ) );
end