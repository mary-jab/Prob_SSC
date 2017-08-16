%
% function CMat = CVX_Lasso(Xp, Q, cst, lambda0, lambda1)
%
% if (nargin < 3)
%     cst = 0;
% end
%
% if (nargin < 4)
%     lambda0 = 0.001;
% end
% if (nargin < 5)
%     lambda1 = 10*lambda0;
% end
% D = size(Xp,1);
% N = size(Xp,2);
% tic
% for i = 1:N
%     y = Xp(:,i);
%     q = Q(:,i); q(i) = [];
%     if i == 1
%         Y = Xp(:,i+1:end);
%     elseif ( (i > 1) && (i < N) )
%         Y = [Xp(:,1:i-1) Xp(:,i+1:N)];
%     else
%         Y = Xp(:,1:N-1);
%     end
%
%     % L1 optimization using CVX
%     if cst == 1
%         cvx_begin quiet;
%         cvx_precision high
%         variable c(N-1,1);
%         minimize( norm(c,1) +lambda0 * norm(Y * c  - y) +lambda1* norm((1-q).*c,2));
%         subject to
%         sum(c) == 1;
%         cvx_end;
%     else
%         cvx_begin quiet;
%         cvx_precision high
%         variable c(N-1,1);
%         minimize( norm(c,1) + lambda0 * norm(Y * c  - y)+lambda1* norm((1-q).*c,2));
%         cvx_end;
%     end
%
%     % place 0's in the diagonals of the coefficient matrix
%     if i == 1
%         CMat(1,1) = 0;
%         CMat(2:N,1) = c(1:N-1);
%     elseif ( (i > 1) && (i < N) )
%         CMat(1:i-1,i) = c(1:i-1);
%         CMat(i,i) = 0;
%         CMat(i+1:N,i) = c(i:N-1);
%     else
%         CMat(1:N-1,N) = c(1:N-1);
%         CMat(N,N) = 0;
%     end
%
% end
% toc
% end

% %% ssc
% function [ZZ] = CVX_Lasso(Y, QMat,cst, lambda0, lambda1)
% tic
% N = (size(Y,2));
% cvx_begin quiet;
% %cvx_precision high
% variable ZZ(N,N);
% expression X(N)
% for i = 1: N,
%     X(i)=(norm (Y * ZZ(:,i)  - Y(:,i),2)) + lambda0* norm(ZZ(:,i),1)+...
%         lambda1*( norm((1-QMat(:,i)).* ZZ(:,i),2)   ) ;
% end
% minimize(sum(X))
% subject to
% diag(ZZ) == 0;
% cvx_end;
% toc
% end

%% vectorized
function [ZZ] = CVX_Lasso(Y, QMat,cst, lambda0, lambda1)
Opt =  'L1Noisy';
tic   
N = (size(Y,2));

if ( strcmp(Opt , 'Lasso') )
    cvx_begin quiet;
    %cvx_precision high
    variable ZZ(N,N);
    expression X(N)
    for i = 1: N,
        X(i)=(norm (Y * ZZ(:,i)  - Y(:,i),2)) + lambda0* norm(ZZ(:,i),1)+...
            lambda1*( norm((1-QMat(:,i)).* ZZ(:,i),2)   ) ;
    end
    minimize(sum(X))
    subject to
    diag(ZZ) == 0;
    cvx_end;
    
elseif ( strcmp(Opt , 'L1Noisy') )
    cvx_begin;
    cvx_precision high
    expression X(N)
    variable ZZ(N,N);
    for i = 1: N,
        X(i)=lambda0* norm(ZZ(:,i),1)+...
            lambda1*( norm((1-QMat(:,i)).* ZZ(:,i),2)   ) ;
    end
    minimize(sum(X))
    subject to
    norm (Y * ZZ  - Y,2) <= lambda0+lambda1;
    diag(ZZ) == 0;
    cvx_end;
    
end




toc
end

% function [ZZ] = CVX_Lasso(Y, QMat,cst, lambda0, lambda1)
% tic
%
% N = (size(Y,2));
% diagIDX = 1:(N+1):N*N;
%
% % Now vectorize everything
% I = speye(N);
% J = kron(I,Y);
% y = Y(:);
%
% cvx_begin quiet;
% %cvx_precision high
% variable ZZ(N*N);
% expressions X
% % for i = 1: N,
%     X=(norm (J * ZZ  - y,2)) + lambda0* norm(ZZ,1);%+...
% %         lambda1*( norm((1-QMat(:,i)).* ZZ(:,i),2)   ) ;
% % end
% minimize(X)
% subject to
% ZZ(diagIDX)==0;
% cvx_end;
%
% ZZ = reshape(ZZ,[N,N]);
% ZZ(abs(ZZ)<.0001)=0;
% toc
% end


%
% function [ZZ] = CVX_Lasso(Y, QMat,cst, lambda0, lambda1)
% tic
% N = (size(Y,2));
% cvx_begin quiet;
% %cvx_precision high
% variables ZZ(N,N) Err zz aa;
%
% minimize(Err + lambda0 * zz + lambda1*aa)
% subject to
% norm(Y*ZZ - Y, 2) <= Err;
% norm(ZZ, 1) <= zz;
% norm((1-QMat).* ZZ, 'fro') <=aa
% diag(ZZ) == 0;
%
% cvx_end;
% toc
% end
