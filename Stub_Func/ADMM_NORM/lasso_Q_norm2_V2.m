%--------------------------------------------------------------------------
% This function takes a DxN matrix of N data points in a D-dimensional
% space and returns a NxN coefficient matrix of the sparse representation
% of each data point in terms of the rest of the points
% Y: DxN data matrix
% affine: if true then enforce the affine constraint
% thr1: stopping threshold for the coefficient error ||Z-C||
% thr2: stopping threshold for the linear system error ||Y-YZ||
% maxIter: maximum number of iterations of ADMM
% C2: NxN sparse coefficient matrix
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------

function C2 = lasso_Q_norm2_V2(Y, lmbd0, lmbd1, Q, preZ, alpha,gamma0 , options, maxIter,thr)

[D,N] = size(Y);
if (nargin < 4)
    Q = ones(N,N);
end
if (nargin < 5)
    preZ = zeros(N,N);
end
if (nargin < 6)
    % default regularizarion parameters
    alpha = 800;
end
if (nargin < 7)
    options.ouliers = 0;
end
if (nargin < 8)
    gamma0 =0.1;
end
if (nargin < 9)
    % default maximum number of iterations of ADMM
    maxIter = 200;
end
if (nargin < 10)
    % default coefficient error threshold to stop ADMM
    % default linear system error threshold to stop ADMM
    thr = 2*10^-4;
end
if (length(alpha) == 1)
    alpha1 = alpha(1);
    alpha2 = alpha(1);
    alpha3 = alpha(1);
elseif (length(alpha) == 2)
    alpha1 = alpha(1);
    alpha2 = alpha(2);
    alpha3 = alpha(2);
elseif (length(alpha) == 3)
    alpha1 = alpha(1);
    alpha2 = alpha(2);
    alpha3 = alpha(3);
end

if (length(thr) == 1)
    thr1 = thr(1);
    thr2 = thr(1);
elseif (length(thr) == 2)
    thr1 = thr(1);
    thr2 = thr(2);
end

% setting penalty parameters for the ADMM
mu1  = alpha1*1/computeLambda_mat(Y);%1/alpha1 * lmbd0;%lmbd0;%
mu11 =  mu1*lmbd1;%1/alpha1 * lmbd1;lmbd1;%
mu2 = alpha2 * 1;


if (~options.ouliers)
    % initialization
    A = inv(mu1*(Y'*Y)+ mu11*(1-Q)*(1-Q)' +mu2*eye(N));
    C1 = preZ;
    Lambda2 = zeros(N,N);
    err1 = 10*thr1; err2 = 10*thr2;
    i = 1;
    % ADMM iterations
    while ( err1(i) > thr1 && i < maxIter )
        % updating Z
        Z = A * (mu1*(Y'*Y)+mu2*(C1-Lambda2/mu2));
        Z = Z - diag(diag(Z));
        % updating C
        C2 = max(0,(abs(Z+Lambda2/mu2) - 1/mu2*ones(N))) .* sign(Z+Lambda2/mu2);
        C2 = C2 - diag(diag(C2));
        % updating Lagrange multipliers
        Lambda2 = Lambda2 + mu2 * (Z - C2);
        % computing errors
        err1(i+1) = errorCoef(Z,C2);
        err2(i+1) = errorLinSys(Y,Z);
        %
        C1 = C2;
        i = i + 1;
    end
    fprintf('err1: %2.4f, err2: %2.4f, iter: %3.0f \n',err1(end),err2(end),i);
else
    
    gamma = alpha3 / norm(Y,1);
    P = [Y eye(D)/gamma];
    Qext = [Q eye(N,D)/gamma ];
    
    A = inv(mu1*(P'*P)+ mu11*(1-Qext)'*(1-Qext) +mu2*eye(N+D));
    C1 =[preZ; zeros(D, N)];
    Lambda1 = zeros(D,N);
    Lambda2 = zeros(N+D,N);
    err1 = 10*thr1; err2 = 10*thr2;
    i = 1;
    % ADMM iterations
    while ( (err1(i) > thr1 || err2(i) > thr2) && i < maxIter )
        % updating Z
        Z = A * (mu1*P'*(Y+Lambda1/mu1)+mu2*(C1-Lambda2/mu2));
        Z(1:N,:) = Z(1:N,:) - diag(diag(Z(1:N,:)));
        % updating C
        C2(N+1:N+D,:)  = max(0,(abs(Z(N+1:N+D,:)+Lambda2(N+1:N+D,:) /mu2) - 1/mu2*ones(D,N))) .* sign(Z(N+1:N+D,:)+Lambda2(N+1:N+D,:) /mu2);
        C2(1:N,:)      = max(0,(abs(Z(1:N,:)    +Lambda2(1:N,:)/mu2)      - 1/mu2*ones(N,N))) .*  sign(Z(1:N,:)+Lambda2(1:N,:)/mu2);
        
        C2(1:N,:) = C2(1:N,:) - diag(diag(C2(1:N,:)));
        % updating Lagrange multipliers
        Lambda1 = Lambda1 + mu1 * (Y - P * Z);
        Lambda2 = Lambda2 + mu2 * (Z - C2);
        % computing errors
        err1(i+1) = errorCoef(Z,C2);
        err2(i+1) = errorLinSys(P,Z);
        %
        C1 = C2;
        i = i + 1;
    end
    fprintf('err1: %2.4f, err2: %2.4f, iter: %3.0f \n',err1(end),err2(end),i);
    C2 =C2(1:N,:);
end
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


